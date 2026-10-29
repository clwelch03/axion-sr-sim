#!/usr/bin/env python3
"""Pull this account's Codex chat stores into the current cluster node.

Python 3.11+; standard library only. Preview by default, --apply to import.
See docs/codex-chat-recovery.md before changing the storage assumptions.
"""

import argparse
import contextlib
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import select
import shlex
import shutil
import socket
import sqlite3
import subprocess
import sys
import tempfile
import time
import tomllib


DATABASES = ("state_5.sqlite", "thread_history_1.sqlite")
TABLES = {
    "threads": (0, "id"),
    "thread_dynamic_tools": (0, "thread_id"),
    "thread_attachments": (0, "thread_id"),
    "thread_spawn_edges": (0, "child_thread_id"),
    "guardian_review_feedback": (0, "thread_id"),
    "thread_read_receipts": (0, "thread_id"),
    "thread_turns": (1, "thread_id"),
    "thread_items": (1, "thread_id"),
    "thread_history_projection_state": (1, "thread_id"),
    "thread_realtime_items": (1, "thread_id"),
}
# Older sources predate these additive tables and contribute no rows for them.
# All shared table definitions must match; newer sources cannot be imported
# into a destination that lacks their tables.
OPTIONAL_TABLES = {"guardian_review_feedback", "thread_read_receipts"}
CONTENT_TABLES = tuple(t for t in TABLES if t not in (
    "threads", "thread_history_projection_state", "thread_read_receipts"))
HISTORY_TABLES = tuple(t for t, (database, _) in TABLES.items() if database == 1)
ORGANIZATION = ("name", "project_id", "thread_section_id", "section_position",
                "section_entered_at_ms", "is_pinned", "archived", "archived_at")
# Recency can change when a node merely opens a chat. Import these timestamps
# along with real changes, but do not schedule a recovery solely for UI activity.
RECENCY = ("updated_at", "updated_at_ms", "recency_at", "recency_at_ms")
UUID_PATTERN = re.compile(r"[0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12}")
DEFAULT_NODES = tuple("n000%d.scs00" % i for i in (1, 2, 3))


class RecoveryError(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise RecoveryError(message)


def connect(path, writable=False):
    require(path.is_file(), "Database missing: %s" % path)
    connection = sqlite3.connect(path.resolve().as_uri() +
                                 ("?mode=rw" if writable else "?mode=ro"),
                                 uri=True, timeout=20)
    connection.row_factory = sqlite3.Row
    return connection


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def snapshot(source, destination, timeout=60):
    """Use SQLite's backup API, including committed WAL data; never copy live DBs."""
    destination.mkdir(mode=0o700, parents=True, exist_ok=False)
    # Keep SQLite's writable scratch files on local storage, not the shared mount.
    with contextlib.ExitStack() as stack:
        temporary = stack.enter_context(tempfile.TemporaryDirectory(prefix="codex-chat-snapshot-"))
        readers = [stack.enter_context(contextlib.closing(connect(source / filename))) for filename in DATABASES]
        versions = [reader.execute("PRAGMA data_version").fetchone()[0] for reader in readers]
        for filename, reader in zip(DATABASES, readers):
            local = Path(temporary) / filename
            with contextlib.closing(sqlite3.connect(local)) as writer:
                deadline = time.monotonic() + timeout

                def progress(status, remaining, total):
                    require(time.monotonic() < deadline, "Snapshot timed out: " + filename)

                reader.backup(writer, pages=256, progress=progress)
                require(writer.execute("PRAGMA quick_check").fetchone()[0] == "ok",
                        "Invalid snapshot: " + filename)
                writer.execute("PRAGMA journal_mode=DELETE")
        require(versions == [reader.execute("PRAGMA data_version").fetchone()[0] for reader in readers],
                "Chat databases changed during snapshot; rerun when writes settle.")
        # Check the captured pair before copying it to potentially slow shared
        # storage. Later source writes do not invalidate these completed backups.
        for filename in DATABASES:
            shutil.copyfile(Path(temporary) / filename, destination / filename)


class Store:
    """An immutable in-memory view of a saved pair of databases."""

    def __init__(self, directory):
        self.directory = directory
        self.schema = {}
        self.columns = {}
        self.rows = {}
        self.primary_keys = {}
        for number, filename in enumerate(DATABASES):
            with contextlib.closing(connect(directory / filename)) as connection:
                all_tables = {r[0] for r in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'")}
                for table in all_tables - TABLES.keys():
                    columns = {r[1] for r in connection.execute('PRAGMA table_info("%s")' % table)}
                    require(not columns.intersection({"thread_id", "child_thread_id", "parent_thread_id"}),
                            "Unsupported thread table: " + table)
                for table, (database, group_key) in TABLES.items():
                    if database != number:
                        continue
                    if table in OPTIONAL_TABLES and table not in all_tables:
                        continue
                    require(table in all_tables, "Missing table: " + table)
                    info = list(connection.execute("PRAGMA table_info(" + table + ")"))
                    self.columns[table] = [r[1] for r in info]
                    self.primary_keys[table] = [r[1] for r in sorted(info, key=lambda r: r[5]) if r[5]]
                    self.schema[table] = connection.execute(
                        "SELECT sql FROM sqlite_master WHERE name=?", (table,)).fetchone()[0]
                    groups = {}
                    for result in connection.execute("SELECT * FROM " + table):
                        row = dict(result)
                        key = tuple(row[k] for k in self.primary_keys[table])
                        groups.setdefault(row[group_key], {})[key] = row
                    self.rows[table] = groups
        self.ids = set().union(*(set(groups) for groups in self.rows.values()))

    def records(self, table, thread_id):
        return self.rows.get(table, {}).get(thread_id, {})

    def thread(self, thread_id):
        return next(iter(self.records("threads", thread_id).values()), None)

    def contains(self, other, thread_id):
        complete = all(all(self.records(table, thread_id).get(key) == row
                       for key, row in other.records(table, thread_id).items())
                   for table in CONTENT_TABLES)
        for key, row in other.records("thread_history_projection_state", thread_id).items():
            ours = self.records("thread_history_projection_state", thread_id).get(key)
            complete = complete and ours is not None and all(
                ours[column] >= row[column] for column in ("next_rollout_byte_offset", "next_rollout_ordinal"))
        return complete


def continuation_boundary(path, base):
    """Check a continuation's byte/ordinal boundary and inventory its old tail.

    A continuation can deliberately branch before the end of its parent. Keep
    that whole parent (and its database rows), and disclose the omitted tail;
    never rewrite the boundary to splice two different histories together.
    """
    end_ordinal, end_byte = base.get("end_ordinal_exclusive"), base.get("end_byte_offset")
    require(type(end_ordinal) is int and type(end_byte) is int
            and end_ordinal >= 0 and end_byte >= 0,
            "Invalid continuation boundary: " + str(path))
    offset, boundary_offset, tail_records, tail_turns = 0, None, 0, set()
    with path.open("rb") as stream:
        for index, line in enumerate(stream):
            record = json.loads(line)
            if index == 0:
                start = (record.get("payload", {}).get("history_base") or {}).get("end_ordinal_exclusive", 0)
                require(type(start) is int and start >= 0, "Invalid ancestor ordinal: " + str(path))
            ordinal = start + index
            if ordinal == end_ordinal:
                boundary_offset = offset
            if ordinal >= end_ordinal:
                tail_records += 1
                payload = record.get("payload", {})
                if record.get("type") == "turn_context" or (
                        record.get("type") == "event_msg" and payload.get("type") == "task_started"):
                    if payload.get("turn_id"):
                        tail_turns.add(payload["turn_id"])
            offset += len(line)
        if offset and ordinal + 1 == end_ordinal:
            boundary_offset = offset
    require(boundary_offset == end_byte,
            "Continuation byte/ordinal boundary does not match ancestor: " + str(path))
    return {"path": str(path.resolve()), "end_ordinal_exclusive": end_ordinal,
            "end_byte_offset": end_byte, "tail_records": tail_records,
            "tail_turn_ids": sorted(tail_turns)}


def continuation_paths(thread_id, selected, paths, home, *, audit_boundaries=True):
    """Verify the selected rollout descends from every known indexed version."""
    require(home is not None, "A Codex home is required to verify continuations")
    inventory = transcript_inventory(home)
    pending = Path(selected)
    hashes, names, boundaries = {}, set(), []
    child_base = None
    while True:
        path = pending.resolve()
        require(path.is_relative_to(home.resolve()), "Rollout outside shared Codex home: " + str(path))
        require(str(path) not in hashes, "Cyclic continuation history: " + thread_id)
        require(path.is_file(), "Missing rollout: " + str(path))
        hashes[str(path)] = digest(path)
        names.add(path.name)
        with path.open() as stream:
            first = json.loads(next(stream))
        require(first.get("type") == "session_meta", "Unsupported rollout header: " + str(path))
        metadata = first["payload"]
        require(metadata.get("id") == thread_id, "Continuation belongs to another chat: " + str(path))
        if child_base is not None and audit_boundaries:
            boundaries.append(continuation_boundary(path, child_base))
        base = metadata.get("history_base")
        if not base:
            break
        parent = base.get("thread_id")
        matches = inventory.get(parent, [])
        require(matches, "Missing ancestor rollout: " + str(parent))
        require(len({digest(match) for match in matches}) == 1, "Conflicting ancestor rollouts: " + str(parent))
        pending = matches[0]
        child_base = base
    require({Path(path).name for path in paths} <= names,
            "Selected continuation does not descend from every known rollout: " + thread_id)
    check_transcripts(hashes)
    return hashes, boundaries


def uncached_history_path(source, thread_id, home):
    """Identify a complete rollout whose pagination cache can be rebuilt.

    Never fall back from a partial cache or a continuation: those need their
    database history. The isolated engine must still read the imported chat.
    """
    row = source.thread(thread_id)
    if home is None or not row or row.get("history_mode") != "paginated":
        return None
    if any(source.records(table, thread_id) for table in (
            "thread_turns", "thread_items", "thread_realtime_items", "thread_history_projection_state")):
        return None
    path = Path(row["rollout_path"]).resolve()
    require(path.is_relative_to(home.resolve()), "Rollout outside shared Codex home: " + str(path))
    expected = digest(path)
    with path.open() as stream:
        first = json.loads(stream.readline())
    require(first.get("type") == "session_meta", "Unsupported rollout header: " + str(path))
    metadata = first["payload"]
    require(metadata.get("id") == thread_id, "Rollout belongs to another chat: " + str(path))
    if metadata.get("history_base") or UUID_PATTERN.findall(path.name)[-1:] != [thread_id]:
        return None
    if thread_id in empty_helper_threads(source, [thread_id]):
        return None
    check_transcripts({str(path): expected})
    return path, expected


def make_plan(stores, excluded, reviewed_continuations=None, home=None):
    """Choose a source containing every known content row, never by timestamp alone."""
    target = stores["local"]
    reviewed_continuations = reviewed_continuations or {}
    for thread_id, label in reviewed_continuations.items():
        require(label in stores and stores[label].thread(thread_id),
                "Reviewed continuation source has no indexed chat: %s=%s" % (thread_id, label))
    for label, store in stores.items():
        mismatches = [table for table in TABLES
                      if store.schema.get(table) != target.schema.get(table)
                      and not (table in OPTIONAL_TABLES and table not in store.schema)]
        require(not mismatches, "Schema mismatch: %s (%s)" % (label, ", ".join(mismatches)))
    # Continuation history can live under an ID present only in a rollout filename.
    excluded = set(excluded)
    for thread_id in list(excluded):
        for store in stores.values():
            row = store.thread(thread_id)
            if row:
                excluded.update(UUID_PATTERN.findall(Path(row["rollout_path"]).name))
    plan, conflicts = [], []
    for thread_id in sorted(set().union(*(store.ids for store in stores.values())) - excluded):
        existing = target.thread(thread_id)
        known = [store for store in stores.values() if thread_id in store.ids]
        # An indexed chat must be sourced from an indexed version.
        indexed = any(store.thread(thread_id) for store in known)
        candidates = [(label, store) for label, store in stores.items()
                      if thread_id in store.ids and (not indexed or store.thread(thread_id))
                      and all(store.contains(other, thread_id) for other in known)]
        if not candidates:
            conflicts.append({"id": thread_id, "reason": "Divergent content; no source contains all versions"})
            continue
        # Choose lineage by ancestry, never by timestamps. A stale node can still
        # index an earlier segment of a chat resumed elsewhere.
        paths = {store.thread(thread_id)["rollout_path"] for store in known if store.thread(thread_id)}
        # Archiving moves a file without changing its identity.
        verified_paths, boundaries = None, []
        if len({Path(path).name for path in paths}) > 1:
            reviewed_source = reviewed_continuations.get(thread_id)
            if reviewed_source is not None:
                candidates = [(label, store) for label, store in candidates if label == reviewed_source]
                if not candidates:
                    conflicts.append({"id": thread_id, "reason": "Reviewed source does not contain every known content row"})
                    continue
            verified, failures = {}, []
            for label, store in candidates:
                try:
                    verified[label] = continuation_paths(
                        thread_id, store.thread(thread_id)["rollout_path"], paths, home,
                        audit_boundaries=reviewed_source is None)
                except (RecoveryError, OSError, ValueError) as error:
                    if reviewed_source is not None:
                        raise
                    failures.append("%s: %s" % (label, error))
            candidates = [(label, store) for label, store in candidates if label in verified]
            if not candidates:
                conflicts.append({"id": thread_id,
                                  "reason": "Different rollout/continuation paths could not be verified; " + "; ".join(failures),
                                  "title": next(store.thread(thread_id).get("title", thread_id)
                                                for store in known if store.thread(thread_id)),
                                  "paths": {label: store.thread(thread_id)["rollout_path"]
                                            for label, store in stores.items() if store.thread(thread_id)}})
                continue
        def rank(candidate):
            label, store = candidate
            row = store.thread(thread_id) or {}
            return (row.get("history_mode") == "paginated",
                    row.get("updated_at_ms") or row.get("updated_at", 0) * 1000, label == "local")
        label, source = max(candidates, key=rank)
        if len({Path(path).name for path in paths}) > 1:
            verified_paths, boundaries = verified[label]
        source_row = source.thread(thread_id)
        content_changed = any(source.records(table, thread_id) != target.records(table, thread_id)
                              for table in CONTENT_TABLES)
        metadata_changed = source_row and existing and any(
            source_row.get(key) != existing.get(key)
            for key in source.columns["threads"] if key not in ORGANIZATION + RECENCY)
        if label == "local" or (not content_changed and not metadata_changed and bool(source_row) == bool(existing)):
            continue
        uncached = uncached_history_path(source, thread_id, home)
        # New project/section memberships need additional tables; fail explicitly.
        if source_row and not existing and (source_row.get("project_id") or source_row.get("thread_section_id")):
            conflicts.append({"id": thread_id, "reason": "New chat belongs to a project or section; manual import required"})
            continue
        plan.append({"id": thread_id, "node": label,
                     "kind": "new" if source_row and not existing else "updated" if indexed else "continuation",
                     "title": (source_row or {}).get("name") or (source_row or {}).get("title", thread_id),
                     "source": (source_row or {}).get("source"),
                     "metadata_only": not content_changed and bool(existing)})
        if verified_paths:
            plan[-1]["verified_continuation_paths"] = verified_paths
            plan[-1]["continuation_boundaries"] = boundaries
        if uncached is not None:
            path, expected = uncached
            plan[-1]["rebuild_history_paths"] = {str(path): expected}
    return plan, conflicts, excluded


def desired_rows(stores, entry, table):
    thread_id = entry["id"]
    if table in HISTORY_TABLES and "recovered_history" in entry:
        return {tuple(row[k] for k in stores["local"].primary_keys[table]): dict(row)
                for row in entry["recovered_history"][table]}
    rows = {key: dict(row) for key, row in stores[entry["node"]].records(table, thread_id).items()}
    existing = stores["local"].thread(thread_id)
    # Read/unread cursors are local UI state, not competing chat histories.
    # Preserve them (including absence) for existing chats; seed new chats from
    # the source unless the destination already has a receipt for that ID.
    if table == "thread_read_receipts":
        local_rows = stores["local"].records(table, thread_id)
        if existing or local_rows:
            return {key: dict(row) for key, row in local_rows.items()}
    if table == "threads" and existing:
        for row in rows.values():
            for column in ORGANIZATION:
                if column != "name" or existing[column] is not None:
                    row[column] = existing[column]
    return rows


def merge(destination, stores, plan):
    """Compare touched records with the snapshot under write locks, then import."""
    original = stores["local"]
    with contextlib.closing(connect(destination / DATABASES[0], writable=True)) as connection:
        connection.execute("ATTACH DATABASE ? AS history", (str(destination / DATABASES[1]),))
        connection.execute("PRAGMA foreign_keys=ON")
        connection.execute("BEGIN IMMEDIATE")
        try:
            for table, (database, group_key) in TABLES.items():
                schema = "history" if database else "main"
                sql = connection.execute("SELECT sql FROM " + schema + ".sqlite_master WHERE name=?", (table,)).fetchone()
                if table not in original.schema:
                    require(sql is None, "Destination schema changed: " + table)
                    continue
                require(sql and sql[0] == original.schema[table], "Destination schema changed: " + table)
                qualified = schema + "." + table
                columns = original.columns[table]
                for entry in plan:
                    thread_id = entry["id"]
                    existing = {tuple(row[k] for k in original.primary_keys[table]): dict(row)
                                for row in connection.execute("SELECT * FROM " + qualified + " WHERE " + group_key + "=?", (thread_id,))}
                    require(existing == original.records(table, thread_id),
                            "Destination changed since snapshot: %s / %s. Rerun." % (thread_id, table))
                    desired = desired_rows(stores, entry, table)
                    if table == "threads" and existing:
                        row = next(iter(desired.values()))
                        connection.execute("UPDATE " + qualified + " SET " + ",".join(c + "=?" for c in columns) + " WHERE id=?",
                                           [row[c] for c in columns] + [thread_id])
                    else:
                        connection.execute("DELETE FROM " + qualified + " WHERE " + group_key + "=?", (thread_id,))
                        connection.executemany("INSERT INTO " + qualified + " VALUES (" + ",".join("?" for _ in columns) + ")",
                                               [[row[c] for c in columns] for row in desired.values()])
            for schema in ("main", "history"):
                require(not connection.execute("PRAGMA " + schema + ".foreign_key_check").fetchall(),
                        "Foreign key check failed: " + schema)
            connection.commit()
        except BaseException:
            connection.rollback()
            raise


def verify_rows(destination, stores, plan):
    actual = Store(destination)
    for entry in plan:
        for table in TABLES:
            require(actual.records(table, entry["id"]) == desired_rows(stores, entry, table),
                    "Post-import mismatch: %s / %s" % (entry["id"], table))
    for filename in DATABASES:
        with contextlib.closing(connect(destination / filename)) as connection:
            require(connection.execute("PRAGMA quick_check").fetchone()[0] == "ok", "Integrity failure: " + filename)


def transcript_inventory(home):
    inventory = {}
    for folder in ("sessions", "archived_sessions"):
        for path in (home / folder).rglob("*.jsonl"):
            ids = UUID_PATTERN.findall(path.name)
            if ids:
                inventory.setdefault(ids[-1], []).append(path)
    return inventory


def collect_transcripts(stores, plan, home):
    inventory = transcript_inventory(home)
    indexed = {}
    for store in stores.values():
        for thread_id in store.rows["threads"]:
            indexed[thread_id] = store.thread(thread_id)
    # Read the paths that will actually be installed, not whichever store was
    # visited last (which can still point to an older continuation).
    indexed.update({tid: stores["local"].thread(tid) for tid in stores["local"].rows["threads"]})
    for entry in plan:
        row = stores[entry["node"]].thread(entry["id"])
        if row:
            indexed[entry["id"]] = row
        check_transcripts(entry.get("verified_continuation_paths", {}))
        check_transcripts(entry.get("rebuild_history_paths", {}))
    verify_ids = {entry["id"] for entry in plan if entry["id"] in indexed}
    for entry in plan:
        if entry["id"] not in indexed:
            owners = {tid for tid, row in indexed.items() if entry["id"] in Path(row["rollout_path"]).name}
            # An intermediate continuation need not appear in the latest path.
            # Its session header retains the owning conversation's ID.
            for path in inventory.get(entry["id"], []):
                with path.open() as stream:
                    first = json.loads(next(stream))
                require(first.get("type") == "session_meta", "Unsupported rollout header: " + str(path))
                owner = first["payload"].get("id")
                if owner in indexed:
                    owners.add(owner)
            require(owners, "Cannot find the owning chat for continuation " + entry["id"])
            verify_ids.update(owners)
    pending = [Path(indexed[tid]["rollout_path"]) for tid in verify_ids]
    pending.extend(path for entry in plan for path in inventory.get(entry["id"], []))
    pending.extend(Path(path) for entry in plan for path in entry.get("verified_continuation_paths", {}))
    hashes = {}
    while pending:
        path = pending.pop().resolve()
        require(path.is_relative_to(home.resolve()), "Rollout outside shared Codex home: " + str(path))
        require(path.is_file(), "Missing rollout: " + str(path))
        if str(path) in hashes:
            continue
        hashes[str(path)] = digest(path)
        with path.open() as stream:
            first = json.loads(next(stream))
        require(first.get("type") == "session_meta", "Unsupported rollout header: " + str(path))
        metadata = first["payload"]
        parents = {metadata.get("forked_from_id"), (metadata.get("history_base") or {}).get("thread_id")}
        for parent in parents - {None}:
            matches = inventory.get(parent, [])
            require(matches, "Missing ancestor rollout: " + parent)
            require(len({digest(match) for match in matches}) == 1, "Conflicting ancestor rollouts: " + parent)
            pending.append(matches[0])
    return sorted(verify_ids), hashes


def check_transcripts(hashes):
    for path, expected in hashes.items():
        require(digest(Path(path)) == expected, "Transcript changed; rerun retrieval: " + path)


def find_codex(explicit):
    if explicit:
        binary = Path(explicit).expanduser().resolve()
    else:
        # Prefer the engine running this session to another installed extension version.
        pid = os.getpid()
        binary = None
        while pid > 1:
            try:
                executable = Path("/proc/%d/exe" % pid).resolve()
                if executable.name == "codex":
                    binary = executable
                    break
                pid = int(next(line.split()[1] for line in Path("/proc/%d/status" % pid).read_text().splitlines()
                               if line.startswith("PPid:")))
            except (OSError, StopIteration):
                break
        if binary is None:
            candidates = list((Path.home() / ".vscode-server/extensions").glob("openai.chatgpt-*/bin/linux-x86_64/codex"))
            binary = max(candidates, key=lambda path: path.stat().st_mtime) if candidates else Path(shutil.which("codex") or "/missing-codex")
    require(binary.is_file() and os.access(binary, os.X_OK), "Codex executable not found; pass --codex-bin.")
    return binary


def empty_helper_threads(store, thread_ids):
    """Recognize guardian placeholders with no turns in either storage format."""
    empty = set()
    for thread_id in thread_ids:
        row = store.thread(thread_id)
        if not row or row.get("source") != '{"subagent":{"other":"guardian"}}':
            continue
        if any(store.records(table, thread_id) for table in
               ("thread_turns", "thread_items", "thread_realtime_items")):
            continue
        with Path(row["rollout_path"]).open() as stream:
            records = [json.loads(line) for line in stream]
        if not records or records[0].get("type") != "session_meta":
            continue
        metadata = records[0]["payload"]
        if metadata.get("id") != thread_id or metadata.get("history_base") or metadata.get("forked_from_id"):
            continue
        if all((record.get("type"), record.get("payload", {}).get("type")) in
               {("response_item", "compaction"), ("event_msg", "thread_settings_applied")}
               for record in records[1:]):
            empty.add(thread_id)
    return empty


def validate_app_server(binary, probe, thread_ids, report_path, resume_ids=()):
    """Read isolated histories; resuming copies starts no model turns."""
    store = Store(probe / "db")
    empty_helpers = empty_helper_threads(store, thread_ids)
    environment = dict(os.environ, CODEX_HOME=str(probe / "home"), CODEX_SQLITE_HOME=str(probe / "db"))
    environment.pop("CODEX_THREAD_ID", None)
    with report_path.with_suffix(".stderr.log").open("w") as errors:
        process = subprocess.Popen([str(binary), "-c", "sqlite_home=" + json.dumps(str(probe / "db")),
                                    "app-server", "--stdio"], env=environment, stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE, stderr=errors, bufsize=0)
        pending = b""
        request_id = 0

        def request(method, params):
            nonlocal pending, request_id
            request_id += 1
            process.stdin.write((json.dumps({"id": request_id, "method": method, "params": params}) + "\n").encode())
            deadline = time.monotonic() + 45
            while time.monotonic() < deadline:
                if b"\n" not in pending:
                    if not select.select([process.stdout], [], [], 1)[0]:
                        continue
                    chunk = os.read(process.stdout.fileno(), 1024 * 1024)
                    require(chunk, "Validation app-server exited; see " + str(errors.name))
                    pending += chunk
                    continue
                line, pending = pending.split(b"\n", 1)
                response = json.loads(line)
                if response.get("id") == request_id:
                    require("error" not in response, "%s: %s" % (method, response.get("error")))
                    return response["result"]
            raise RecoveryError("Validation timed out: " + method)

        results = []
        try:
            request("initialize", {"clientInfo": {"name": "chat_recovery_check", "version": "1"},
                                   "capabilities": {"experimentalApi": True}})
            process.stdin.write(b'{"method":"initialized"}\n')
            for thread_id in thread_ids:
                # This engine builds missing projections on resume. Legacy
                # readers also require a loaded session for list_turns. These
                # writes are confined to disposable copies of the databases
                # and rollouts; never resume a live/source chat here.
                if thread_id in resume_ids or store.thread(thread_id).get("history_mode") == "legacy":
                    request("thread/resume", {"threadId": thread_id})
                # This installed engine supports includeTurns for both storage modes.
                # If a later engine removes it, stop rather than accepting an empty read.
                thread = request("thread/read", {"threadId": thread_id, "includeTurns": True})["thread"]
                require(thread.get("id") == thread_id and "turns" in thread,
                        "Validation did not return the requested chat history: " + thread_id)
                turns = thread["turns"]
                require(turns or thread_id in empty_helpers, "Validation returned no turns: " + thread_id)
                results.append({"id": thread_id, "turn_ids": [turn["id"] for turn in turns], "turns": len(turns)})
                if not turns:
                    results[-1]["empty_helper"] = True
                print("Validated %s: %d turns" % (thread_id, len(turns)), flush=True)
            save_json(report_path, {"binary": str(binary), "threads": results})
        finally:
            process.stdin.close()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            process.stdout.close()


def prepare_probe(probe, root, stores, plan, home, timeout):
    snapshot(root / "local", probe / "db", timeout=timeout)
    merge(probe / "db", stores, plan)
    verify_rows(probe / "db", stores, plan)
    shutil.copytree(root / "transcripts", probe / "home")
    with contextlib.closing(connect(probe / "db" / DATABASES[0], writable=True)) as connection:
        for row in connection.execute("SELECT id,rollout_path FROM threads").fetchall():
            path = Path(row["rollout_path"]).resolve()
            require(path.is_relative_to(home), "Rollout outside Codex home: " + str(path))
            connection.execute("UPDATE threads SET rollout_path=? WHERE id=?",
                               (str(probe / "home" / path.relative_to(home)), row["id"]))
        connection.commit()


def recovered_history(store, thread_id, original, resumed):
    """Accept only rows derived from the original bytes, excluding resume settings."""
    require(resumed.startswith(original), "History rebuild rewrote the copied rollout: " + thread_id)
    for line in resumed[len(original):].splitlines():
        record = json.loads(line)
        require(record.get("type") == "event_msg"
                and record.get("payload", {}).get("type") == "thread_settings_applied",
                "History rebuild appended conversation content: " + thread_id)
    size, count = len(original), len(original.splitlines())
    rows = {table: [dict(row) for row in store.records(table, thread_id).values()] for table in HISTORY_TABLES}
    require(rows["thread_turns"] and len(rows["thread_history_projection_state"]) == 1,
            "History rebuild produced no complete projection: " + thread_id)
    require(not rows["thread_realtime_items"], "Unexpected realtime history during rebuild: " + thread_id)
    for table in ("thread_turns", "thread_items"):
        for row in rows[table]:
            for column in ("rollout_byte_offset", "rollout_end_byte_offset"):
                require(row.get(column) is None or 0 <= row[column] <= size,
                        "Rebuilt history extends beyond original rollout: " + thread_id)
            for column in ("rollout_ordinal", "rollout_end_ordinal", "updated_at_ordinal"):
                require(row.get(column) is None or 0 <= row[column] < count,
                        "Rebuilt history extends beyond original ordinals: " + thread_id)
    projection = rows["thread_history_projection_state"][0]
    require(projection["next_rollout_byte_offset"] == len(resumed)
            and projection["next_rollout_ordinal"] == len(resumed.splitlines()),
            "History rebuild did not consume the full rollout: " + thread_id)
    projection.update(next_rollout_byte_offset=size, next_rollout_ordinal=count)
    return rows


def rebuild_histories(binary, root, stores, plan, home, timeout):
    entries = [entry for entry in plan if entry.get("rebuild_history_paths")]
    if not entries:
        return
    with tempfile.TemporaryDirectory(prefix="codex-chat-rebuild-") as temporary:
        probe = Path(temporary)
        prepare_probe(probe, root, stores, plan, home, timeout)
        ids = [entry["id"] for entry in entries]
        validate_app_server(binary, probe, ids, root / "rebuild-validation.json", resume_ids=ids)
        rebuilt = Store(probe / "db")
        require(rebuilt.schema == stores["local"].schema, "History rebuild changed the database schema")
        for entry in entries:
            path = Path(next(iter(entry["rebuild_history_paths"]))).relative_to(home)
            entry["recovered_history"] = recovered_history(
                rebuilt, entry["id"], (root / "transcripts" / path).read_bytes(),
                (probe / "home" / path).read_bytes())


@contextlib.contextmanager
def writer_locks(home, thread_ids):
    """Exclude conversation writers; SQLite transactions serialize DB writes."""
    with contextlib.ExitStack() as stack:
        lock_directory = home / "thread-writer-locks"
        require(lock_directory.is_dir(), "Missing Codex writer-lock directory")
        # .sqlite-maintenance.lock belongs to a long-lived background worker,
        # including while idle. It is not a lock for ordinary database writes.
        # Taking it would prevent recovery from any running Codex session.
        paths = [lock_directory / ".coordination.lock"]
        paths.extend(lock_directory / (thread_id + ".lock") for thread_id in sorted(thread_ids))
        for path in paths:
            handle = stack.enter_context(path.open("a+b"))
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise RecoveryError("Chat writer or recovery coordination is active; retry later: " + path.name) from error
        yield


def canonical_node(value):
    if re.fullmatch(r"n\d{4}", value):
        value += ".scs00"
    require(bool(re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9.-]*", value)), "Invalid node name: " + value)
    return value


def remote_snapshot(node, command, timeout):
    """Keep remote failures actionable without silently omitting a source."""
    try:
        subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", "-o", "LogLevel=ERROR",
                        node, shlex.join(command)], check=True, timeout=timeout)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
        if isinstance(error, subprocess.TimeoutExpired):
            reason = "Remote snapshot timed out after %s seconds on %s" % (timeout, node)
        elif error.returncode == 255:
            reason = "SSH failed for %s; check the SSH diagnostic above" % node
        else:
            reason = "Remote snapshot failed on %s (exit %d)" % (node, error.returncode)
        raise RecoveryError(
            reason + ". No chats imported; completed backups are retained. "
            "Retry when the source is available, explicitly select reachable sources with --nodes, "
            "or use --from-snapshot LABEL=DIR for saved backups. "
            "Omitting a node leaves its unsaved history unchecked.") from error


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true", help="Import after snapshots and validation; otherwise preview only")
    parser.add_argument("--nodes", nargs="+", help="Source nodes (default n0001..n0003, excluding current node)")
    parser.add_argument("--from-snapshot", action="append", default=[], metavar="LABEL=DIR", help="Use saved source DBs instead of SSH (repeatable)")
    parser.add_argument("--codex-home", type=Path, default=Path(os.environ.get("CODEX_HOME", str(Path.home() / ".codex"))))
    parser.add_argument("--sqlite-home", type=Path, help="Destination DB directory; default env/config sqlite_home")
    parser.add_argument("--backup-root", type=Path, default=Path(__file__).resolve().parents[1] / "recovery/chat-retrieval")
    parser.add_argument("--codex-bin", help="Override automatic discovery of the installed Codex engine")
    parser.add_argument("--remote-python", default="python3.12")
    parser.add_argument("--snapshot-timeout", type=int, default=60, metavar="SECONDS",
                        help="Time limit per database backup (default 60; increase for slow shared storage)")
    parser.add_argument("--exclude-thread", action="append", default=[], help="Also exclude CODEX_THREAD_ID automatically")
    parser.add_argument("--reviewed-continuation", action="append", default=[], metavar="THREAD_ID=SOURCE",
                        help="Use a manually reviewed continuation source; still checks ancestry and all content rows")
    parser.add_argument("--snapshot-to", type=Path, help=argparse.SUPPRESS)
    return parser.parse_args()


def run(args):
    os.umask(0o077)
    require(args.snapshot_timeout > 0, "--snapshot-timeout must be positive.")
    remote_timeout = 2 * args.snapshot_timeout + 30
    home = args.codex_home.expanduser().resolve()
    config_path = home / "config.toml"
    config = tomllib.loads(config_path.read_text()) if config_path.exists() else {}
    configured = os.environ.get("CODEX_SQLITE_HOME") or config.get("sqlite_home")
    require(args.sqlite_home or configured, "sqlite_home is not configured; pass --sqlite-home explicitly.")
    live = Path(args.sqlite_home or configured).expanduser().resolve()
    if args.snapshot_to:
        snapshot(live, args.snapshot_to, timeout=args.snapshot_timeout)
        return
    require(not (args.nodes and args.from_snapshot), "Choose --nodes or --from-snapshot.")
    args.backup_root.mkdir(mode=0o700, parents=True, exist_ok=True)
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    root = Path(tempfile.mkdtemp(prefix=stamp + "-", dir=args.backup_root.resolve()))
    print("Recovery directory: " + str(root), flush=True)
    snapshot(live, root / "local", timeout=args.snapshot_timeout)
    stores = {"local": Store(root / "local")}
    source_commands = {}
    if args.from_snapshot:
        sources = [value.split("=", 1) for value in args.from_snapshot]
        require(all(len(value) == 2 for value in sources), "Expected --from-snapshot LABEL=DIR")
        for label, path in sources:
            require(re.fullmatch(r"[a-zA-Z0-9_-]+", label) and label not in stores, "Invalid or duplicate snapshot label")
            snapshot(Path(path).expanduser().resolve(), root / label, timeout=args.snapshot_timeout)
            stores[label] = Store(root / label)
    else:
        local_node = canonical_node(socket.gethostname())
        for node in dict.fromkeys(canonical_node(value) for value in (args.nodes or DEFAULT_NODES)):
            if node == local_node:
                continue
            label = node.replace(".", "_")
            command = [args.remote_python, str(Path(__file__).resolve()), "--snapshot-to", str(root / label),
                       "--sqlite-home", str(live), "--codex-home", str(home),
                       "--snapshot-timeout", str(args.snapshot_timeout)]
            print("Snapshotting " + node, flush=True)
            remote_snapshot(node, command, remote_timeout)
            source_commands[label] = (node, command)
            stores[label] = Store(root / label)
    excluded = set(args.exclude_thread)
    if os.environ.get("CODEX_THREAD_ID"):
        excluded.add(os.environ["CODEX_THREAD_ID"])
    reviewed = {}
    for value in args.reviewed_continuation:
        parts = value.split("=", 1)
        require(len(parts) == 2 and all(parts), "Expected --reviewed-continuation THREAD_ID=SOURCE")
        thread_id, label = parts
        require(thread_id not in reviewed, "Duplicate reviewed continuation: " + thread_id)
        reviewed[thread_id] = label
    plan, conflicts, excluded = make_plan(stores, excluded, reviewed, home)
    report = {"destination": str(live), "host": socket.gethostname(), "sources": list(stores)[1:],
              "excluded": sorted(excluded), "threads": plan, "conflicts": conflicts}
    save_json(root / "plan.json", report)
    for entry in plan:
        print("Planned %s: %s [%s]" % (entry["kind"], entry["title"][:100], entry["node"]), flush=True)
        if entry.get("rebuild_history_paths"):
            print("  Missing pagination cache; rebuilding and validating from an isolated rollout copy.", flush=True)
        for boundary in entry.get("continuation_boundaries", []):
            if boundary["tail_records"]:
                print("  Continuation branches before %d older records (%d turn starts); "
                      "retaining the full ancestor and all stored history. See plan.json."
                      % (boundary["tail_records"], len(boundary["tail_turn_ids"])), flush=True)
    for conflict in conflicts:
        print("Conflict: %s [%s]: %s" % (conflict.get("title", conflict["id"]),
                                        conflict["id"], conflict["reason"]), flush=True)
    require(not conflicts, "%d conflicts; no import performed. See %s" % (len(conflicts), root / "plan.json"))
    if not plan:
        save_json(root / "receipt.json", {"status": "up_to_date", "changed": 0})
        print("Already up to date; no live databases changed.")
        return
    verify_ids, hashes = collect_transcripts(stores, plan, home)
    require(not (set(verify_ids) & excluded), "Import depends on the current/excluded chat; retry in a fresh chat.")
    save_json(root / "transcripts.json", hashes)
    for path in hashes:
        relative = Path(path).relative_to(home)
        destination = root / "transcripts" / relative
        destination.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
        shutil.copyfile(path, destination)
        require(digest(destination) == hashes[path], "Transcript changed while backing up: " + path)
    binary = find_codex(args.codex_bin)
    rebuild_histories(binary, root, stores, plan, home, args.snapshot_timeout)
    save_json(root / "plan.json", report)
    with tempfile.TemporaryDirectory(prefix="codex-chat-probe-") as temporary:
        probe = Path(temporary)
        prepare_probe(probe, root, stores, plan, home, args.snapshot_timeout)
        validate_app_server(binary, probe, verify_ids, root / "validation.json")
    check_transcripts(hashes)
    if not args.apply:
        print("Preview validated. Rerun with --apply for a fresh audit and import.")
        return
    # This also locks continuation ancestor IDs to protect history dependencies.
    lock_ids = set(verify_ids) | {entry["id"] for entry in plan}
    for path in hashes:
        lock_ids.update(UUID_PATTERN.findall(Path(path).name))
    require(not (lock_ids & excluded), "Import depends on an excluded conversation; retry in a fresh chat.")
    with writer_locks(home, lock_ids):
        check_transcripts(hashes)
        # Catch source updates made between its first DB snapshot and hashing the
        # shared rollout files. Fresh snapshots are read-only on the source node.
        for label, (node, command) in source_commands.items():
            check_directory = root / (label + "-recheck")
            command = list(command)
            command[command.index("--snapshot-to") + 1] = str(check_directory)
            remote_snapshot(node, command, remote_timeout)
            fresh = Store(check_directory)
            require(fresh.schema == stores[label].schema, "Source schema changed: " + label)
            for thread_id in lock_ids:
                for table in TABLES:
                    require(fresh.records(table, thread_id) == stores[label].records(table, thread_id),
                            "Source changed during retrieval: %s / %s. Rerun." % (label, thread_id))
        check_transcripts(hashes)
        snapshot(live, root / "before-apply", timeout=args.snapshot_timeout)
        save_json(root / "apply-started.json", {"status": "applying", "backup": str(root / "before-apply")})
        merge(live, stores, plan)
        verify_rows(live, stores, plan)
        check_transcripts(hashes)
    counts = {kind: sum(entry["kind"] == kind for entry in plan) for kind in ("new", "updated", "continuation")}
    save_json(root / "receipt.json", {"status": "applied", "counts": counts, "backup": str(root / "before-apply")})
    print("Import verified: " + json.dumps(counts))
    print("Reload the VS Code window if the chats do not appear immediately.")


def main():
    try:
        run(parse_args())
    except (RecoveryError, OSError, ValueError, sqlite3.Error, subprocess.SubprocessError) as error:
        print("Recovery stopped: " + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
