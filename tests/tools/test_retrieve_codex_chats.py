"""Data preservation checks for the cluster chat recovery utility."""

import fcntl
import importlib.util
from pathlib import Path
import sqlite3
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
