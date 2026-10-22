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
    "thread_turns": (1, "thread_id"),
    "thread_items": (1, "thread_id"),
    "thread_history_projection_state": (1, "thread_id"),
    "thread_realtime_items": (1, "thread_id"),
}
CONTENT_TABLES = tuple(t for t in TABLES if t not in (
    "threads", "thread_history_projection_state"))
ORGANIZATION = ("name", "project_id", "thread_section_id", "section_position",
                "section_entered_at_ms", "is_pinned", "archived", "archived_at")
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


def snapshot(source, destination):
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
                deadline = time.monotonic() + 60

                def progress(status, remaining, total):
                    require(time.monotonic() < deadline, "Snapshot timed out: " + filename)

                reader.backup(writer, pages=256, progress=progress)
                require(writer.execute("PRAGMA quick_check").fetchone()[0] == "ok",
                        "Invalid snapshot: " + filename)
                writer.execute("PRAGMA journal_mode=DELETE")
            shutil.copyfile(local, destination / filename)
        require(versions == [reader.execute("PRAGMA data_version").fetchone()[0] for reader in readers],
                "Chat databases changed during snapshot; rerun when writes settle.")


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
        return self.rows[table].get(thread_id, {})

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


def make_plan(stores, excluded):
    """Choose a source containing every known content row, never by timestamp alone."""
    target = stores["local"]
    for label, store in stores.items():
        require(store.schema == target.schema, "Schema mismatch: " + label)
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
        # Do not silently change continuation lineage, even with newer timestamps.
        paths = {store.thread(thread_id)["rollout_path"] for store in known if store.thread(thread_id)}
        # Archiving moves a file without changing its identity.
        if len({Path(path).name for path in paths}) > 1:
            conflicts.append({"id": thread_id, "reason": "Different rollout/continuation paths require manual review"})
            continue
        def rank(candidate):
            label, store = candidate
            row = store.thread(thread_id) or {}
            return (row.get("history_mode") == "paginated",
                    row.get("updated_at_ms") or row.get("updated_at", 0) * 1000, label == "local")
        label, source = max(candidates, key=rank)
        source_row = source.thread(thread_id)
        content_changed = any(source.records(table, thread_id) != target.records(table, thread_id)
                              for table in CONTENT_TABLES)
        metadata_changed = source_row and existing and any(
            source_row.get(key) != existing.get(key)
            for key in source.columns["threads"] if key not in ORGANIZATION)
        if label == "local" or (not content_changed and not metadata_changed and bool(source_row) == bool(existing)):
            continue
        # New project/section memberships need additional tables; fail explicitly.
        if source_row and not existing and (source_row.get("project_id") or source_row.get("thread_section_id")):
            conflicts.append({"id": thread_id, "reason": "New chat belongs to a project or section; manual import required"})
            continue
        plan.append({"id": thread_id, "node": label,
                     "kind": "new" if source_row and not existing else "updated" if indexed else "continuation",
                     "title": (source_row or {}).get("name") or (source_row or {}).get("title", thread_id),
                     "source": (source_row or {}).get("source"),
                     "metadata_only": not content_changed and bool(existing)})
    return plan, conflicts, excluded


def desired_rows(stores, entry, table):
    thread_id = entry["id"]
    rows = {key: dict(row) for key, row in stores[entry["node"]].records(table, thread_id).items()}
    existing = stores["local"].thread(thread_id)
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
    verify_ids = {entry["id"] for entry in plan if entry["id"] in indexed}
    for entry in plan:
        if entry["id"] not in indexed:
            owners = {tid for tid, row in indexed.items() if entry["id"] in Path(row["rollout_path"]).name}
            require(owners, "Cannot find the owning chat for continuation " + entry["id"])
            verify_ids.update(owners)
    pending = [Path(indexed[tid]["rollout_path"]) for tid in verify_ids]
    pending.extend(path for entry in plan for path in inventory.get(entry["id"], []))
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


def validate_app_server(binary, probe, thread_ids, report_path):
    """Read full history in an isolated home. No turns, auth, or source writes."""
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
                # This installed engine supports includeTurns for both storage modes.
                # If a later engine removes it, stop rather than accepting an empty read.
                thread = request("thread/read", {"threadId": thread_id, "includeTurns": True})["thread"]
                turns = thread.get("turns", [])
                require(turns, "Validation returned no turns: " + thread_id)
                results.append({"id": thread_id, "turn_ids": [turn["id"] for turn in turns], "turns": len(turns)})
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


@contextlib.contextmanager
def writer_locks(home, database_home, thread_ids):
    with contextlib.ExitStack() as stack:
        lock_directory = home / "thread-writer-locks"
        require(lock_directory.is_dir(), "Missing Codex writer-lock directory")
        paths = [database_home / ".sqlite-maintenance.lock", lock_directory / ".coordination.lock"]
        paths.extend(lock_directory / (thread_id + ".lock") for thread_id in sorted(thread_ids))
        for path in paths:
            handle = stack.enter_context(path.open("a+b"))
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise RecoveryError("Chat or maintenance is active; retry later: " + path.name) from error
        yield


def canonical_node(value):
    if re.fullmatch(r"n\d{4}", value):
        value += ".scs00"
    require(bool(re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9.-]*", value)), "Invalid node name: " + value)
    return value


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
    parser.add_argument("--exclude-thread", action="append", default=[], help="Also exclude CODEX_THREAD_ID automatically")
    parser.add_argument("--snapshot-to", type=Path, help=argparse.SUPPRESS)
    return parser.parse_args()


def run(args):
    os.umask(0o077)
    home = args.codex_home.expanduser().resolve()
    config_path = home / "config.toml"
    config = tomllib.loads(config_path.read_text()) if config_path.exists() else {}
    configured = os.environ.get("CODEX_SQLITE_HOME") or config.get("sqlite_home")
    require(args.sqlite_home or configured, "sqlite_home is not configured; pass --sqlite-home explicitly.")
    live = Path(args.sqlite_home or configured).expanduser().resolve()
    if args.snapshot_to:
        snapshot(live, args.snapshot_to)
        return
    require(not (args.nodes and args.from_snapshot), "Choose --nodes or --from-snapshot.")
    args.backup_root.mkdir(mode=0o700, parents=True, exist_ok=True)
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    root = Path(tempfile.mkdtemp(prefix=stamp + "-", dir=args.backup_root.resolve()))
    print("Recovery directory: " + str(root), flush=True)
    snapshot(live, root / "local")
    stores = {"local": Store(root / "local")}
    source_commands = {}
    if args.from_snapshot:
        sources = [value.split("=", 1) for value in args.from_snapshot]
        require(all(len(value) == 2 for value in sources), "Expected --from-snapshot LABEL=DIR")
        for label, path in sources:
            require(re.fullmatch(r"[a-zA-Z0-9_-]+", label) and label not in stores, "Invalid or duplicate snapshot label")
            snapshot(Path(path).expanduser().resolve(), root / label)
            stores[label] = Store(root / label)
    else:
        local_node = canonical_node(socket.gethostname())
        for node in dict.fromkeys(canonical_node(value) for value in (args.nodes or DEFAULT_NODES)):
            if node == local_node:
                continue
            label = node.replace(".", "_")
            command = [args.remote_python, str(Path(__file__).resolve()), "--snapshot-to", str(root / label),
                       "--sqlite-home", str(live), "--codex-home", str(home)]
            print("Snapshotting " + node, flush=True)
            subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", "-o", "LogLevel=ERROR",
                            node, shlex.join(command)], check=True, timeout=150)
            source_commands[label] = (node, command)
            stores[label] = Store(root / label)
    excluded = set(args.exclude_thread)
    if os.environ.get("CODEX_THREAD_ID"):
        excluded.add(os.environ["CODEX_THREAD_ID"])
    plan, conflicts, excluded = make_plan(stores, excluded)
    report = {"destination": str(live), "host": socket.gethostname(), "sources": list(stores)[1:],
              "excluded": sorted(excluded), "threads": plan, "conflicts": conflicts}
    save_json(root / "plan.json", report)
    for entry in plan:
        print("%s: %s [%s]" % (entry["kind"], entry["title"][:100], entry["node"]), flush=True)
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
    with tempfile.TemporaryDirectory(prefix="codex-chat-probe-") as temporary:
        probe = Path(temporary)
        snapshot(root / "local", probe / "db")
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
    with writer_locks(home, live, lock_ids):
        check_transcripts(hashes)
        # Catch source updates made between its first DB snapshot and hashing the
        # shared rollout files. Fresh snapshots are read-only on the source node.
        for label, (node, command) in source_commands.items():
            check_directory = root / (label + "-recheck")
            command = list(command)
            command[command.index("--snapshot-to") + 1] = str(check_directory)
            subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", "-o", "LogLevel=ERROR",
                            node, shlex.join(command)], check=True, timeout=150)
            fresh = Store(check_directory)
            require(fresh.schema == stores[label].schema, "Source schema changed: " + label)
            for thread_id in lock_ids:
                for table in TABLES:
                    require(fresh.records(table, thread_id) == stores[label].records(table, thread_id),
                            "Source changed during retrieval: %s / %s. Rerun." % (label, thread_id))
        check_transcripts(hashes)
        snapshot(live, root / "before-apply")
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
