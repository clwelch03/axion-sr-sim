"""Data preservation checks for the cluster chat recovery utility."""

import fcntl
import importlib.util
import io
import json
import os
from pathlib import Path
import shlex
import sqlite3
import subprocess
import tempfile
import unittest
from unittest import mock


spec = importlib.util.spec_from_file_location(
    "chat_recovery", Path(__file__).resolve().parents[2] / "scripts/retrieve_codex_chats.py")
recovery = importlib.util.module_from_spec(spec)
spec.loader.exec_module(recovery)

STATE_SCHEMA = """
CREATE TABLE threads (
    id TEXT PRIMARY KEY, rollout_path TEXT, title TEXT, name TEXT, source TEXT,
    updated_at INTEGER, updated_at_ms INTEGER, history_mode TEXT,
    project_id TEXT, thread_section_id TEXT, section_position INTEGER,
    section_entered_at_ms INTEGER, is_pinned INTEGER DEFAULT 0,
    archived INTEGER DEFAULT 0, archived_at INTEGER);
CREATE TABLE thread_dynamic_tools (thread_id TEXT REFERENCES threads(id), position INTEGER,
    name TEXT, PRIMARY KEY(thread_id,position));
CREATE TABLE thread_attachments (id TEXT PRIMARY KEY, thread_id TEXT REFERENCES threads(id), payload TEXT);
CREATE TABLE thread_spawn_edges (parent_thread_id TEXT, child_thread_id TEXT PRIMARY KEY, status TEXT);
CREATE TABLE guardian_review_feedback (
    id TEXT PRIMARY KEY NOT NULL,
    thread_id TEXT NOT NULL REFERENCES threads(id) ON DELETE CASCADE,
    record BLOB NOT NULL
);
CREATE TABLE thread_read_receipts (
    thread_id TEXT PRIMARY KEY NOT NULL,
    first_unread_turn TEXT,
    revision TEXT NOT NULL,
    last_published_turn TEXT
);
"""
HISTORY_SCHEMA = """
CREATE TABLE thread_turns (thread_id TEXT, turn_id TEXT, status TEXT, PRIMARY KEY(thread_id,turn_id));
CREATE TABLE thread_items (thread_id TEXT, turn_id TEXT, item_id TEXT, item_json TEXT,
    PRIMARY KEY(thread_id,turn_id,item_id));
CREATE TABLE thread_realtime_items (thread_id TEXT, item_id TEXT, item_json TEXT, PRIMARY KEY(thread_id,item_id));
CREATE TABLE thread_history_projection_state (thread_id TEXT PRIMARY KEY,
    next_rollout_byte_offset INTEGER, next_rollout_ordinal INTEGER);
"""


def create(directory):
    directory.mkdir()
    for filename, schema in zip(recovery.DATABASES, (STATE_SCHEMA, HISTORY_SCHEMA)):
        with sqlite3.connect(directory / filename) as connection:
            connection.executescript(schema)


def add_chat(directory, thread_id="chat", timestamp=1, name="Title"):
    with sqlite3.connect(directory / recovery.DATABASES[0]) as connection:
        connection.execute("INSERT INTO threads(id,rollout_path,title,name,source,updated_at,history_mode) VALUES (?,?,?,?,?,?,?)",
                           (thread_id, "/home/test/.codex/sessions/" + thread_id + ".jsonl", "Initial prompt", name, "vscode", timestamp, "paginated"))


def add_turn(directory, thread_id="chat", turn="turn", content="message"):
    with sqlite3.connect(directory / recovery.DATABASES[1]) as connection:
        connection.execute("INSERT INTO thread_turns VALUES (?,?,?)", (thread_id, turn, "completed"))
        connection.execute("INSERT INTO thread_items VALUES (?,?,?,?)", (thread_id, turn, "item", content))


class RetrievalTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.local, self.remote = self.root / "local", self.root / "remote"
        create(self.local)
        create(self.remote)

    def stores(self):
        return {"local": recovery.Store(self.local), "remote": recovery.Store(self.remote)}

    def test_import_and_repeat_preserve_organization_and_other_chat(self):
        add_chat(self.local, name="My name")
        add_chat(self.local, "unrelated")
        add_chat(self.remote, timestamp=2, name="Remote name")
        add_chat(self.remote, "new", timestamp=2)
        add_turn(self.local)
        add_turn(self.remote)
        add_turn(self.remote, turn="second")
        add_turn(self.remote, thread_id="new")
        with sqlite3.connect(self.local / recovery.DATABASES[0]) as connection:
            connection.execute("UPDATE threads SET archived=1,archived_at=5,is_pinned=1 WHERE id='chat'")
        stores = self.stores()
        plan, conflicts, _ = recovery.make_plan(stores, set())
        self.assertEqual(conflicts, [])
        # A different chat can receive a turn after planning without being overwritten.
        add_turn(self.local, thread_id="unrelated", content="ongoing work")
        recovery.merge(self.local, stores, plan)
        recovery.verify_rows(self.local, stores, plan)
        actual = recovery.Store(self.local)
        self.assertEqual(actual.thread("chat")["name"], "My name")
        self.assertEqual(actual.thread("chat")["archived"], 1)
        self.assertEqual(actual.thread("chat")["is_pinned"], 1)
        self.assertEqual(len(actual.records("thread_turns", "chat")), 2)
        self.assertEqual(len(actual.records("thread_items", "unrelated")), 1)
        self.assertEqual(recovery.make_plan(self.stores(), set())[:2], ([], []))

    def test_divergent_message_does_not_choose_newest_timestamp(self):
        add_chat(self.local)
        add_chat(self.remote, timestamp=100)
        add_turn(self.local, content="local response")
        add_turn(self.remote, content="different response")
        plan, conflicts, _ = recovery.make_plan(self.stores(), set())
        self.assertEqual(plan, [])
        self.assertEqual(len(conflicts), 1)

    def test_split_new_turns_need_manual_merge(self):
        add_chat(self.local)
        add_chat(self.remote)
        add_turn(self.local, turn="local-only")
        add_turn(self.remote, turn="remote-only")
        self.assertEqual(len(recovery.make_plan(self.stores(), set())[1]), 1)

    def test_destination_change_rolls_back_the_entire_merge(self):
        add_chat(self.local, timestamp=1)
        add_chat(self.remote, timestamp=2)
        add_turn(self.remote)
        stores = self.stores()
        plan, _, _ = recovery.make_plan(stores, set())
        add_turn(self.local, content="arrived after snapshot")
        with self.assertRaisesRegex(recovery.RecoveryError, "Destination changed"):
            recovery.merge(self.local, stores, plan)
        self.assertEqual(recovery.Store(self.local).thread("chat")["updated_at"], 1)
        self.assertEqual(next(iter(recovery.Store(self.local).records("thread_items", "chat").values()))["item_json"], "arrived after snapshot")

    def test_continuation_history_without_index_row_is_imported(self):
        add_turn(self.remote, thread_id="unindexed-continuation")
        stores = self.stores()
        plan, conflicts, _ = recovery.make_plan(stores, set())
        self.assertEqual(conflicts, [])
        self.assertEqual(plan[0]["kind"], "continuation")
        recovery.merge(self.local, stores, plan)
        recovery.verify_rows(self.local, stores, plan)
        self.assertIsNone(recovery.Store(self.local).thread("unindexed-continuation"))
        self.assertEqual(len(recovery.Store(self.local).records("thread_turns", "unindexed-continuation")), 1)

    def test_current_chat_and_its_continuation_are_excluded(self):
        parent = "00000000-0000-0000-0000-000000000001"
        branch = "00000000-0000-0000-0000-000000000002"
        add_chat(self.remote, parent)
        with sqlite3.connect(self.remote / recovery.DATABASES[0]) as connection:
            connection.execute("UPDATE threads SET rollout_path=?", ("/home/test/.codex/sessions/rollout-" + parent + "_" + branch + ".jsonl",))
        add_turn(self.remote, thread_id=branch)
        plan, conflicts, excluded = recovery.make_plan(self.stores(), {parent})
        self.assertEqual((plan, conflicts), ([], []))
        self.assertIn(branch, excluded)

    def test_schema_mismatch_is_rejected(self):
        with sqlite3.connect(self.remote / recovery.DATABASES[0]) as connection:
            connection.execute("ALTER TABLE threads ADD COLUMN new_field TEXT")
        with self.assertRaisesRegex(recovery.RecoveryError, "Schema mismatch"):
            recovery.make_plan(self.stores(), set())

    def test_unhandled_thread_table_is_rejected(self):
        with sqlite3.connect(self.remote / recovery.DATABASES[0]) as connection:
            connection.execute("CREATE TABLE thread_notes (thread_id TEXT, note TEXT)")
        with self.assertRaisesRegex(recovery.RecoveryError, "Unsupported thread table"):
            recovery.Store(self.remote)

    def test_feedback_and_new_chat_receipts_are_imported_exactly(self):
        add_chat(self.local)
        for directory in (self.local, self.remote):
            add_chat(directory, "unrelated")
            with sqlite3.connect(directory / recovery.DATABASES[0]) as connection:
                connection.execute("INSERT INTO guardian_review_feedback VALUES (?,?,?)",
                                   ("unrelated-review", "unrelated", b"unchanged"))
        add_chat(self.remote)
        add_chat(self.remote, "new")
        for thread_id in ("chat", "new"):
            with sqlite3.connect(self.remote / recovery.DATABASES[0]) as connection:
                connection.execute("INSERT INTO guardian_review_feedback VALUES (?,?,?)",
                                   (thread_id + "-review", thread_id, b"\x00\xfffeedback"))
                connection.execute("INSERT INTO thread_read_receipts VALUES (?,?,?,?)",
                                   (thread_id, "unread", "revision", "published"))
        stores = self.stores()
        plan, conflicts, _ = recovery.make_plan(stores, set())
        self.assertEqual(conflicts, [])
        self.assertEqual({entry["id"] for entry in plan}, {"chat", "new"})
        recovery.merge(self.local, stores, plan)
        recovery.verify_rows(self.local, stores, plan)
        actual = recovery.Store(self.local)
        for thread_id in ("chat", "new", "unrelated"):
            self.assertEqual(actual.records("guardian_review_feedback", thread_id),
                             stores["remote"].records("guardian_review_feedback", thread_id))
        self.assertEqual(actual.records("thread_read_receipts", "chat"), {})
        self.assertEqual(actual.records("thread_read_receipts", "new"),
                         stores["remote"].records("thread_read_receipts", "new"))
        self.assertEqual(recovery.make_plan(self.stores(), set())[:2], ([], []))

    def test_divergent_feedback_does_not_choose_newest_timestamp(self):
        for directory, timestamp, record in ((self.local, 1, b"local"), (self.remote, 100, b"remote")):
            add_chat(directory, timestamp=timestamp)
            with sqlite3.connect(directory / recovery.DATABASES[0]) as connection:
                connection.execute("INSERT INTO guardian_review_feedback VALUES ('review','chat',?)", (record,))
        plan, conflicts, _ = recovery.make_plan(self.stores(), set())
        self.assertEqual(plan, [])
        self.assertEqual(len(conflicts), 1)

    def test_different_read_receipts_preserve_destination_state(self):
        for directory, revision in ((self.local, "local"), (self.remote, "remote")):
            add_chat(directory)
            with sqlite3.connect(directory / recovery.DATABASES[0]) as connection:
                connection.execute("INSERT INTO thread_read_receipts VALUES ('chat',NULL,?,NULL)", (revision,))
        self.assertEqual(recovery.make_plan(self.stores(), set())[:2], ([], []))
        add_turn(self.remote)
        stores = self.stores()
        plan, conflicts, _ = recovery.make_plan(stores, set())
        self.assertEqual(conflicts, [])
        self.assertEqual(len(plan), 1)
        recovery.merge(self.local, stores, plan)
        recovery.verify_rows(self.local, stores, plan)
        self.assertEqual(recovery.Store(self.local).records("thread_read_receipts", "chat"),
                         stores["local"].records("thread_read_receipts", "chat"))

    def test_changed_feedback_or_receipts_roll_back_import(self):
        add_chat(self.local)
        add_chat(self.remote, timestamp=100)
        add_turn(self.remote)
        for table, statement in (
            ("guardian_review_feedback", "INSERT INTO guardian_review_feedback VALUES ('review','chat',X'00FF')"),
            ("thread_read_receipts", "INSERT INTO thread_read_receipts VALUES ('chat',NULL,'revision',NULL)"),
        ):
            with self.subTest(table=table):
                stores = self.stores()
                plan, conflicts, _ = recovery.make_plan(stores, set())
                self.assertEqual(conflicts, [])
                with sqlite3.connect(self.local / recovery.DATABASES[0]) as connection:
                    connection.execute(statement)
                with self.assertRaisesRegex(recovery.RecoveryError, "Destination changed.*" + table):
                    recovery.merge(self.local, stores, plan)
                actual = recovery.Store(self.local)
                self.assertEqual(actual.thread("chat")["updated_at"], 1)
                self.assertEqual(actual.records("thread_items", "chat"), {})
                self.assertTrue(actual.records(table, "chat"))
                with sqlite3.connect(self.local / recovery.DATABASES[0]) as connection:
                    connection.execute("DELETE FROM " + table)

    def test_older_matching_stores_without_optional_tables_can_merge(self):
        for directory in (self.local, self.remote):
            with sqlite3.connect(directory / recovery.DATABASES[0]) as connection:
                for table in recovery.OPTIONAL_TABLES:
                    connection.execute("DROP TABLE " + table)
        add_chat(self.remote)
        add_turn(self.remote)
        stores = self.stores()
        plan, conflicts, _ = recovery.make_plan(stores, set())
        self.assertEqual(conflicts, [])
        recovery.merge(self.local, stores, plan)
        recovery.verify_rows(self.local, stores, plan)
        self.assertIsNotNone(recovery.Store(self.local).thread("chat"))

    def test_older_source_imports_into_newer_destination(self):
        add_chat(self.local)
        add_chat(self.remote, timestamp=2)
        add_chat(self.remote, "new")
        add_turn(self.remote)
        add_turn(self.remote, thread_id="new")
        with sqlite3.connect(self.local / recovery.DATABASES[0]) as connection:
            connection.execute("INSERT INTO thread_read_receipts VALUES ('chat',NULL,'local',NULL)")
        with sqlite3.connect(self.remote / recovery.DATABASES[0]) as connection:
            for table in recovery.OPTIONAL_TABLES:
                connection.execute("DROP TABLE " + table)
        stores = self.stores()
        plan, conflicts, _ = recovery.make_plan(stores, set())
        self.assertEqual(conflicts, [])
        self.assertEqual({entry["id"] for entry in plan}, {"chat", "new"})
        recovery.merge(self.local, stores, plan)
        recovery.verify_rows(self.local, stores, plan)
        actual = recovery.Store(self.local)
        self.assertEqual(actual.thread("chat")["updated_at"], 2)
        self.assertEqual(actual.records("thread_read_receipts", "chat"),
                         stores["local"].records("thread_read_receipts", "chat"))
        self.assertEqual(actual.records("thread_read_receipts", "new"), {})
        self.assertEqual(actual.schema, stores["local"].schema)
        self.assertEqual(recovery.Store(self.remote).schema, stores["remote"].schema)
        self.assertEqual(recovery.make_plan(self.stores(), set())[:2], ([], []))

    def test_older_source_cannot_erase_destination_feedback(self):
        add_chat(self.local)
        add_chat(self.remote, timestamp=100)
        with sqlite3.connect(self.local / recovery.DATABASES[0]) as connection:
            connection.execute("INSERT INTO guardian_review_feedback VALUES ('review','chat',X'00FF')")
        with sqlite3.connect(self.remote / recovery.DATABASES[0]) as connection:
            connection.execute("DROP TABLE guardian_review_feedback")
        self.assertEqual(recovery.make_plan(self.stores(), set())[:2], ([], []))
        add_turn(self.remote)
        plan, conflicts, _ = recovery.make_plan(self.stores(), set())
        self.assertEqual(plan, [])
        self.assertEqual(len(conflicts), 1)

    def test_newer_source_cannot_import_into_older_destination(self):
        with sqlite3.connect(self.local / recovery.DATABASES[0]) as connection:
            connection.execute("DROP TABLE guardian_review_feedback")
        with self.assertRaisesRegex(recovery.RecoveryError, r"Schema mismatch: remote \(guardian_review_feedback\)"):
            recovery.make_plan(self.stores(), set())

    def test_changed_optional_table_definition_is_rejected(self):
        with sqlite3.connect(self.remote / recovery.DATABASES[0]) as connection:
            connection.execute("ALTER TABLE guardian_review_feedback ADD COLUMN extra TEXT")
        with self.assertRaisesRegex(recovery.RecoveryError, r"Schema mismatch: remote \(guardian_review_feedback\)"):
            recovery.make_plan(self.stores(), set())

    def test_optional_table_added_after_snapshot_is_rejected(self):
        add_chat(self.remote)
        original_schema = recovery.Store(self.local).schema["guardian_review_feedback"]
        for directory in (self.local, self.remote):
            with sqlite3.connect(directory / recovery.DATABASES[0]) as connection:
                connection.execute("DROP TABLE guardian_review_feedback")
        stores = self.stores()
        plan, _, _ = recovery.make_plan(stores, set())
        with sqlite3.connect(self.local / recovery.DATABASES[0]) as connection:
            connection.execute(original_schema)
        with self.assertRaisesRegex(recovery.RecoveryError, "Destination schema changed: guardian_review_feedback"):
            recovery.merge(self.local, stores, plan)
        self.assertIsNone(recovery.Store(self.local).thread("chat"))

    def test_changed_lineage_and_new_project_are_reported(self):
        add_chat(self.local)
        add_chat(self.remote)
        with sqlite3.connect(self.remote / recovery.DATABASES[0]) as connection:
            connection.execute("UPDATE threads SET rollout_path='/home/test/.codex/sessions/other.jsonl'")
        self.assertIn("continuation paths", recovery.make_plan(self.stores(), set())[1][0]["reason"])
        add_chat(self.remote, "project-chat")
        with sqlite3.connect(self.remote / recovery.DATABASES[0]) as connection:
            connection.execute("UPDATE threads SET project_id='project' WHERE id='project-chat'")
        self.assertEqual(len(recovery.make_plan(self.stores(), set())[1]), 2)

    def continuation_fixture(self):
        owner, middle, latest = ["00000000-0000-0000-0000-%012d" % i for i in (1, 2, 3)]
        home = self.root / "home"
        sessions = home / "sessions"
        sessions.mkdir(parents=True)
        paths = []
        for number, segment in enumerate((owner, middle, latest)):
            name = owner if number == 0 else owner + "_" + segment
            path = sessions / ("rollout-" + name + ".jsonl")
            metadata = {"id": owner, "history_mode": "paginated"}
            if number:
                metadata["history_base"] = {"thread_id": (owner, middle)[number - 1],
                                            "end_ordinal_exclusive": number, "end_byte_offset": 100}
            path.write_text(json.dumps({"type": "session_meta", "payload": metadata}) + "\n")
            paths.append(path)
        for directory, path, timestamp in ((self.local, paths[0], 100), (self.remote, paths[-1], 2)):
            add_chat(directory, owner, timestamp=timestamp)
            with sqlite3.connect(directory / recovery.DATABASES[0]) as connection:
                connection.execute("UPDATE threads SET rollout_path=? WHERE id=?", (str(path), owner))
            add_turn(directory, thread_id=owner)
        for segment in (middle, latest):
            add_turn(self.remote, thread_id=segment)
        return home, (owner, middle, latest), paths

    def test_reviewed_continuation_preserves_all_segments_and_repeat(self):
        home, (owner, middle, latest), paths = self.continuation_fixture()
        stores = self.stores()
        self.assertEqual(len(recovery.make_plan(stores, set())[1]), 1)
        plan, conflicts, _ = recovery.make_plan(stores, set(), {owner: "remote"}, home)
        self.assertEqual(conflicts, [])
        self.assertEqual({entry["id"] for entry in plan}, {owner, middle, latest})
        verify_ids, hashes = recovery.collect_transcripts(stores, plan, home)
        self.assertEqual(verify_ids, [owner])
        self.assertEqual(set(hashes), {str(path) for path in paths})
        recovery.merge(self.local, stores, plan)
        recovery.verify_rows(self.local, stores, plan)
        self.assertEqual(recovery.Store(self.local).thread(owner)["rollout_path"], str(paths[-1]))
        self.assertEqual(recovery.make_plan(self.stores(), set(), {owner: "remote"}, home)[:2], ([], []))

    def test_reviewed_continuation_cannot_choose_older_path_or_skip_content(self):
        home, (owner, _, _), _paths = self.continuation_fixture()
        with self.assertRaisesRegex(recovery.RecoveryError, "does not descend"):
            recovery.make_plan(self.stores(), set(), {owner: "local"}, home)
        add_turn(self.local, thread_id=owner, turn="local-only")
        plan, conflicts, _ = recovery.make_plan(self.stores(), set(), {owner: "remote"}, home)
        self.assertTrue(conflicts)
        self.assertNotIn(owner, {entry["id"] for entry in plan})

    def test_reviewed_continuation_rejects_missing_wrong_and_cyclic_ancestors(self):
        home, (owner, _, latest), paths = self.continuation_fixture()
        for metadata, error in (
            ({"id": owner, "history_base": {"thread_id": "missing"}}, "Missing ancestor"),
            ({"id": "another-chat"}, "another chat"),
            ({"id": owner, "history_base": {"thread_id": latest}}, "Cyclic"),
        ):
            with self.subTest(error=error):
                paths[-1].write_text(json.dumps({"type": "session_meta", "payload": metadata}) + "\n")
                with self.assertRaisesRegex(recovery.RecoveryError, error):
                    recovery.make_plan(self.stores(), set(), {owner: "remote"}, home)

    def test_reviewed_continuation_detects_transcript_change_after_review(self):
        home, (owner, _, _), paths = self.continuation_fixture()
        stores = self.stores()
        plan, _, _ = recovery.make_plan(stores, set(), {owner: "remote"}, home)
        with paths[0].open("a") as stream:
            stream.write('{}\n')
        with self.assertRaisesRegex(recovery.RecoveryError, "Transcript changed"):
            recovery.collect_transcripts(stores, plan, home)

    def test_transcript_collection_uses_selected_source_even_if_listed_first(self):
        home, (owner, _, _), paths = self.continuation_fixture()
        stores = self.stores()
        stores = {"remote": stores["remote"], "local": stores["local"]}
        plan, _, _ = recovery.make_plan(stores, set(), {owner: "remote"}, home)
        verify_ids, hashes = recovery.collect_transcripts(stores, plan, home)
        self.assertEqual(verify_ids, [owner])
        self.assertIn(str(paths[-1]), hashes)

    def test_empty_helper_requires_no_database_or_transcript_turns(self):
        home, (owner, _, _), paths = self.continuation_fixture()
        source = '{"subagent":{"other":"guardian"}}'
        with sqlite3.connect(self.local / recovery.DATABASES[0]) as connection:
            connection.execute("UPDATE threads SET source=?", (source,))
        # A metadata-only file cannot hide already stored history rows.
        self.assertEqual(recovery.empty_helper_threads(recovery.Store(self.local), [owner]), set())
        with sqlite3.connect(self.local / recovery.DATABASES[1]) as connection:
            connection.execute("DELETE FROM thread_turns")
            connection.execute("DELETE FROM thread_items")
        header = paths[0].read_text()
        paths[0].write_text(header + json.dumps({"type": "response_item", "payload": {"type": "compaction"}}) + "\n")
        self.assertEqual(recovery.empty_helper_threads(recovery.Store(self.local), [owner]), {owner})
        for record in (
            {"type": "event_msg", "payload": {"type": "task_started"}},
            {"type": "response_item", "payload": {"type": "message", "role": "user"}},
            {"type": "response_item", "payload": {"type": "function_call"}},
        ):
            with self.subTest(record=record):
                paths[0].write_text(header + json.dumps(record) + "\n")
                self.assertEqual(recovery.empty_helper_threads(recovery.Store(self.local), [owner]), set())
        paths[0].write_text(header)
        with sqlite3.connect(self.local / recovery.DATABASES[0]) as connection:
            connection.execute("UPDATE threads SET source='vscode'")
        self.assertEqual(recovery.empty_helper_threads(recovery.Store(self.local), [owner]), set())

    def test_newer_metadata_cannot_rewind_projection_state(self):
        for directory, timestamp, offset in ((self.local, 1, 100), (self.remote, 20, 50)):
            add_chat(directory, timestamp=timestamp)
            add_turn(directory)
            with sqlite3.connect(directory / recovery.DATABASES[1]) as connection:
                connection.execute("INSERT INTO thread_history_projection_state VALUES ('chat',?,?)", (offset, offset))
        self.assertEqual(recovery.make_plan(self.stores(), set())[:2], ([], []))

    def test_snapshot_contains_committed_wal_rows(self):
        with sqlite3.connect(self.remote / recovery.DATABASES[0]) as writer:
            writer.execute("PRAGMA journal_mode=WAL")
            add_chat(self.remote)
            recovery.snapshot(self.remote, self.root / "snapshot")
        self.assertIsNotNone(recovery.Store(self.root / "snapshot").thread("chat"))

    def test_snapshot_timeout_can_be_extended(self):
        add_chat(self.remote)
        with mock.patch.object(recovery.time, "monotonic", side_effect=range(0, 100, 2)):
            with self.assertRaisesRegex(recovery.RecoveryError, "Snapshot timed out"):
                recovery.snapshot(self.remote, self.root / "timed-out", timeout=1)
        with mock.patch.object(recovery.time, "monotonic", side_effect=range(0, 100, 2)):
            recovery.snapshot(self.remote, self.root / "extended", timeout=10)
        self.assertIsNotNone(recovery.Store(self.root / "extended").thread("chat"))

    def test_source_writes_during_archive_copy_do_not_change_captured_pair(self):
        add_chat(self.remote)
        add_turn(self.remote)
        copyfile = recovery.shutil.copyfile

        def copy_with_later_write(source, destination):
            result = copyfile(source, destination)
            if source.name == recovery.DATABASES[0]:
                add_chat(self.remote, "later")
                add_turn(self.remote, thread_id="later")
            return result

        with mock.patch.object(recovery.shutil, "copyfile", side_effect=copy_with_later_write):
            recovery.snapshot(self.remote, self.root / "snapshot")
        saved = recovery.Store(self.root / "snapshot")
        self.assertEqual(saved.ids, {"chat"})
        self.assertEqual(len(saved.records("thread_turns", "chat")), 1)
        self.assertIn("later", recovery.Store(self.remote).ids)

    def test_source_writes_during_capture_still_reject_snapshot(self):
        add_chat(self.remote)
        connect = recovery.connect

        def connect_with_concurrent_write(path):
            reader = connect(path)
            if path.name != recovery.DATABASES[0]:
                return reader
            wrapped = mock.Mock(wraps=reader)

            def backup_with_write(*args, **kwargs):
                reader.backup(*args, **kwargs)
                add_turn(self.remote, turn="during-snapshot")

            wrapped.backup.side_effect = backup_with_write
            return wrapped

        with mock.patch.object(recovery, "connect", side_effect=connect_with_concurrent_write):
            with self.assertRaisesRegex(recovery.RecoveryError, "changed during snapshot"):
                recovery.snapshot(self.remote, self.root / "snapshot")
        self.assertFalse((self.root / "snapshot" / recovery.DATABASES[0]).exists())

    def test_busy_coordination_lock_stops_import(self):
        home = self.root / "home"
        (home / "thread-writer-locks").mkdir(parents=True)
        with recovery.writer_locks(home, {"chat"}):
            with self.assertRaisesRegex(recovery.RecoveryError, "active"):
                with recovery.writer_locks(home, {"chat"}):
                    self.fail("Concurrent writer was allowed")

    def test_remote_failures_preserve_backups_and_never_apply(self):
        home = self.root / "home"
        (home / "sessions").mkdir(parents=True)
        (home / "thread-writer-locks").mkdir()
        rollout = home / "sessions/chat.jsonl"
        rollout.write_text(json.dumps({"type": "session_meta", "payload": {"id": "chat"}}) + "\n")
        add_chat(self.remote)
        add_turn(self.remote)
        with sqlite3.connect(self.remote / recovery.DATABASES[0]) as connection:
            connection.execute("UPDATE threads SET rollout_path=?", (str(rollout),))
        backup_root = self.root / "backups"
        argv = ["retrieve_codex_chats.py", "--apply", "--nodes", "n0002", "n0003",
                "--codex-home", str(home), "--sqlite-home", str(self.local),
                "--backup-root", str(backup_root)]
        with mock.patch.object(recovery.sys, "argv", argv):
            args = recovery.parse_args()
        original_umask = os.umask(0o077)
        self.addCleanup(os.umask, original_umask)
        for recheck in (False, True):
            for error, message in (
                (subprocess.CalledProcessError(255, "ssh"), "SSH failed for n0003.scs00"),
                (subprocess.CalledProcessError(1, "ssh"), "Remote snapshot failed on n0003.scs00"),
                (subprocess.TimeoutExpired("ssh", 150), "Remote snapshot timed out after 150 seconds on n0003.scs00"),
            ):
                with self.subTest(recheck=recheck, error=message):
                    def ssh(command, **kwargs):
                        remote_command = shlex.split(command[-1])
                        destination = Path(remote_command[remote_command.index("--snapshot-to") + 1])
                        if command[-2] == "n0003.scs00" and destination.name.endswith("-recheck") == recheck:
                            raise error
                        recovery.snapshot(self.remote, destination)
                        return subprocess.CompletedProcess(command, 0)

                    previous = set(backup_root.iterdir()) if backup_root.exists() else set()
                    with (mock.patch.object(recovery.socket, "gethostname", return_value="n0001.scs00"),
                          mock.patch.object(recovery.subprocess, "run", side_effect=ssh),
                          mock.patch.object(recovery, "find_codex", return_value=Path("/unused")),
                          mock.patch.object(recovery, "validate_app_server") as validate,
                          mock.patch.object(recovery, "merge", wraps=recovery.merge) as merge,
                          mock.patch.object(recovery.sys, "stdout", io.StringIO())):
                        with self.assertRaisesRegex(recovery.RecoveryError, message) as caught:
                            recovery.run(args)
                    self.assertIn("No chats imported", str(caught.exception))
                    self.assertIn("--nodes", str(caught.exception))
                    self.assertIn("--from-snapshot", str(caught.exception))
                    run, = set(backup_root.iterdir()) - previous
                    self.assertEqual(recovery.Store(run / "local").ids, set())
                    self.assertIsNotNone(recovery.Store(run / "n0002_scs00").thread("chat"))
                    self.assertEqual(recovery.Store(self.local).ids, set())
                    self.assertFalse((run / "apply-started.json").exists())
                    self.assertFalse((run / "receipt.json").exists())
                    self.assertEqual(validate.call_count, int(recheck))
                    self.assertEqual(merge.call_count, int(recheck))
                    for call in merge.call_args_list:
                        self.assertNotEqual(call.args[0], self.local)
                    with recovery.writer_locks(home, {"chat"}):
                        pass

    def test_busy_thread_lock_stops_import_and_releases_coordination(self):
        home = self.root / "home"
        locks = home / "thread-writer-locks"
        locks.mkdir(parents=True)
        with (locks / "chat.lock").open("a+b") as writer:
            fcntl.flock(writer, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with self.assertRaisesRegex(recovery.RecoveryError, r"active.*chat\.lock"):
                with recovery.writer_locks(home, {"ancestor", "chat"}):
                    self.fail("Active chat was not protected")
            # A failed acquisition must release coordination and earlier locks.
            with recovery.writer_locks(home, {"ancestor"}):
                pass
        with recovery.writer_locks(home, {"chat"}):
            pass

    def test_idle_maintenance_owner_does_not_block_import(self):
        home = self.root / "home"
        locks = home / "thread-writer-locks"
        locks.mkdir(parents=True)
        add_chat(self.local)
        add_chat(self.remote, timestamp=2)
        add_turn(self.remote)
        stores = self.stores()
        plan, conflicts, _ = recovery.make_plan(stores, set())
        self.assertEqual(conflicts, [])
        maintenance_path = self.local / ".sqlite-maintenance.lock"
        with maintenance_path.open("a+b") as maintenance, (locks / "unrelated.lock").open("a+b") as unrelated:
            fcntl.flock(maintenance, fcntl.LOCK_EX | fcntl.LOCK_NB)
            fcntl.flock(unrelated, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with recovery.writer_locks(home, {"chat"}):
                recovery.merge(self.local, stores, plan)
                recovery.verify_rows(self.local, stores, plan)
            # Recovery must not unlink or replace the backend's lock file.
            with maintenance_path.open("rb") as check:
                with self.assertRaises(BlockingIOError):
                    fcntl.flock(check, fcntl.LOCK_EX | fcntl.LOCK_NB)
        self.assertEqual(len(recovery.Store(self.local).records("thread_turns", "chat")), 1)


if __name__ == "__main__":
    unittest.main()
