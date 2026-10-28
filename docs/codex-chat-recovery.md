# Retrieving Codex chats after switching nodes

Use [scripts/retrieve_codex_chats.py](../scripts/retrieve_codex_chats.py) near the
**beginning of a new session**, before resuming a conversation from another node.
Pulling from the new session works even when the previous session ended abruptly.
An end-of-session push would depend on that session staying alive long enough to
run it, and would write to nodes whose chat stores might already be in use.

This is a local recovery utility for this account's Lawrencium setup. It does not
install startup hooks or require a shutdown step. Run it when switching nodes or
when chats or recent messages are missing. It imports into the node on which the
command runs; source chat databases are opened read-only.

## Normal procedure

From this repository, using Python 3.12:

```bash
# Preview: take backups, compare histories, and validate proposed imports.
python3.12 scripts/retrieve_codex_chats.py

# Retrieve: repeat the audit, validate, and import into the current node.
python3.12 scripts/retrieve_codex_chats.py --apply
```

The default sources are `n0001.scs00`, `n0002.scs00`, and `n0003.scs00`, with the
current node excluded automatically. A separate preview is optional: `--apply`
performs all the same checks before writing. Successful imports print
`Import verified`; repeated runs with no new changes print `Already up to date`.
Reload the VS Code window if its chat list has not refreshed.

To choose sources explicitly, for example while working on n0002:

```bash
python3.12 scripts/retrieve_codex_chats.py --nodes n0001 n0003 --apply
```

Short names are expanded to `.scs00`; bare `ssh n0001` did not resolve during the
October 21 recovery. If n0000 has relevant chats, include it explicitly with
`--nodes n0000 n0001 n0002 n0003`.

For future assistant sessions: use this utility instead of adapting the dated
one-off scripts under `recovery/`. A request to retrieve chats authorizes the
import; preview mode is an inspection option, not an extra user-approval step.
If the execution sandbox blocks SSH or writes to Codex's shared writer-lock
directory, request the command's normal sandbox escalation. Do not disable SSH
host checking or change filesystem permissions to work around it.

## Why both files and databases matter

The configuration inspected on October 21, 2026 uses:

| Location | Role |
|---|---|
| `~/.codex/sessions/` and `~/.codex/archived_sessions/` | Shared JSONL rollout files |
| `/tmp/clwelch/codex-sqlite/state_5.sqlite` | Node-local chat index and related records |
| `/tmp/clwelch/codex-sqlite/thread_history_1.sqlite` | Node-local turns, items, and history projection state |
| `~/.codex/thread-writer-locks/` | Shared coordination and conversation writer locks |

Seeing a JSONL file on every node does not mean every node can display the whole
chat. Paginated history also lives in SQLite. Continuations can have history rows
whose IDs have no independent entry in `threads`; copying indexed chats alone
can lose their later messages. Forked chats also depend on ancestor rollouts.

The destination database directory comes from `--sqlite-home`, then
`CODEX_SQLITE_HOME`, then `sqlite_home` in the selected Codex home's `config.toml`.
Codex home defaults to `CODEX_HOME` or `~/.codex`. The script requires an explicit
database location if that setting is absent, rather than guessing from old home
directory databases.

The script assumes the same absolute Codex home, database directory, repository,
and backup paths are available on all selected nodes. The home, repository, and
backups are shared; the SQLite directory is local to each node. It requires
Python 3.11 or newer locally, `python3.12` remotely by default, working batch SSH,
and the installed Codex executable. It uses only the Python standard library.
`--remote-python` and `--codex-bin` override executable discovery when needed.

## What the script does

1. Saves dated SQLite backups from the current node and each requested source.
   It uses SQLite's backup API so committed WAL contents are included. It checks
   database integrity and stops if either database changes during the pair of
   snapshots. It checks that pair on local storage before copying the completed
   backups to shared storage, so later chat activity does not invalidate them.
   It never copies a live `.sqlite` file with `cp` or `rsync`.
2. Compares chat records, turns, items, tools, attachments, spawn edges, guardian
   review feedback, and continuation history. A source must contain every known
   content row with identical values at shared primary keys. A newer timestamp
   alone cannot win over different messages. Projection offsets cannot move backward.
3. Records the proposed changes and any conflicts. It preserves existing
   destination names, pins, archive state, organization, and read/unread receipts.
   A new chat keeps the source's archive state and receipts. Receipt differences
   alone do not trigger imports or history conflicts. Review/helper threads are
   included in counts.
   New project/section memberships and changed continuation paths require manual
   review rather than an incomplete automatic import. A reviewed continuation can
   be selected explicitly as described below; it must descend from every known
   rollout and retain every known database content row.
4. Backs up affected rollouts and their ancestors with SHA-256 hashes. It merges
   into temporary database copies and verifies the exact imported rows. An
   isolated Codex app-server reads every affected chat, including the owners of
   continuation-only history. Validation does not start model turns. A zero-turn
   response is accepted only for a guardian placeholder whose database has no
   history rows and whose rollout contains only its header, initial compaction,
   and settings events, with no ancestor history. Other empty reads stop recovery.
5. With `--apply`, acquires coordination and conversation writer locks, checks the
   transcript hashes again, and takes fresh source snapshots to check that the
   affected records have not changed. It then backs up the destination again.
6. Imports only the selected records. Under SQLite write transactions it checks
   that those destination records still equal the original snapshot. Handled
   errors roll back the transaction. It verifies the result and writes a receipt.

The running `CODEX_THREAD_ID` and its continuation IDs are excluded automatically.
Use `--exclude-thread THREAD_ID` for additional exclusions, or when launching
from a terminal without that environment variable. A busy writer lock stops the
import. The script does not kill sessions or overwrite their rollout files.
The backend's `.sqlite-maintenance.lock` is held by a long-lived background
worker, even while idle; recovery does not acquire or remove it. SQLite write
transactions serialize database changes, while the coordination and per-chat
locks protect the conversations being imported.

The validation API is described in the
[official OpenAI app-server documentation](https://learn.chatgpt.com/docs/app-server).
The installed engine used for the October 21 recovery supports
`thread/read` with `includeTurns: true` for both history modes. Current
documentation deprecates that full-history option for paginated chats in favor
of `thread/turns/list` and `thread/items/list`. If a future engine removes it,
update the validator to page through the full history; do not skip validation.
The SQLite merge itself is based on the locally inspected schema, not a promised
public synchronization API. Shared table definitions must match. Older sources
may omit `guardian_review_feedback` and `thread_read_receipts`; they contribute
no rows for those tables. Existing destination receipts are still preserved,
and local feedback participates in the same content-preservation checks as other
chat records. Other schema mismatches and unknown thread tables stop recovery.

## Backups and receipts

Every run prints its private recovery directory under
`recovery/chat-retrieval/`. Directory names contain a UTC timestamp and a unique
suffix. The directory retains:

| Artifact | Meaning |
|---|---|
| `local/`, source-node directories | Initial database snapshots |
| `plan.json` | Proposed changes, source labels, excluded IDs, and conflicts |
| `transcripts/`, `transcripts.json` | Affected rollouts, ancestors, and hashes |
| `validation.json`, `validation.stderr.log` | Isolated app-server results and diagnostics |
| source directories ending in `-recheck/` | Source snapshots taken just before applying |
| `before-apply/` | Destination backup immediately before the import |
| `apply-started.json` | Marker that the script reached the write phase |
| `receipt.json` | Verified completion, or an already-up-to-date result |

A preview does not write live chat records. It does create recovery artifacts.
Backups contain private conversation data and are excluded from Git. Scratch
space is more durable than a node's `/tmp`, but remains subject to cluster
retention policies. The script never automatically deletes recovery artifacts.

## When something stops the procedure

- **Snapshot timed out:** slow shared storage can exceed the default 60-second
  limit per database backup. Retry with `--snapshot-timeout 300`; the same limit
  is passed to source nodes and their SSH timeout is increased accordingly.
  Backups, integrity checks, and change detection still run.
- **Unsupported thread table: `guardian_review_feedback` or `thread_read_receipts`:**
  older versions of this utility did not handle these tables. Use the corrected
  script and rerun the same command. This failure happens during snapshot
  inspection, before any live import. Feedback is preserved as chat content;
  existing destination read/unread state is retained. Older sources without
  these tables can import into a newer destination when shared table definitions
  match. Importing newer tables into a destination that lacks them still stops.
  Other unknown thread tables still stop recovery for inspection.
- **`.sqlite-maintenance.lock` reported as active:** older versions of this
  utility incorrectly required the backend's background-worker lock. Use the
  corrected script and rerun the same command. That failure occurs before live
  imports; the `new:` and `updated:` lines describe the plan, not applied changes.
  Do not delete the lock file or stop unrelated chats to work around it.
- **Node unavailable or missing databases:** the command fails rather than
  silently treating that node as empty. Retry later, or explicitly select the
  reachable nodes and report that recovery covers only those sources. For
  example, on n0001 with n0003 unavailable, use
  `python3.12 scripts/retrieve_codex_chats.py --nodes n0002 --apply` to check only
  n0002. A connection timeout on SSH port 22 happens before importing chats;
  increasing `--snapshot-timeout` or closing a Codex chat does not fix that
  connection failure. Completed snapshots remain in the printed recovery
  directory. Saved snapshots from an unavailable node can also be used with
  `--from-snapshot` as described below; messages written after them remain
  unchecked. A source failure during the final recheck also stops before any
  live import, so a preview is not sufficient reason to bypass that recheck.
- **Active chat, changed transcript, or database writes during snapshot:** let
  the affected work finish and rerun. Prefer a fresh chat for the recovery so the
  conversation being recovered is not itself excluded or busy.
- **Divergent records:** inspect `plan.json` and saved databases. Do not choose a
  winner by timestamp or use `INSERT OR REPLACE` indiscriminately. Separate
  histories may need preservation as separate conversations.
- **Different rollout/continuation paths:** the terminal prints the conflicting
  title and ID; `plan.json` also lists each source's path. Review the shared
  rollout headers and their `history_base` links, including the content around
  continuation boundaries. After confirming the intended chain, select that
  specific chat and source with
  `--reviewed-continuation THREAD_ID=SOURCE` (repeatable), using the source label
  from `plan.json`, such as `n0003_scs00`. Add this to the normal command with
  `--apply`. The utility still requires the selected source to contain all known
  content rows and verifies that its chain includes every older rollout for the
  same conversation. Missing, cyclic, unrelated, or divergent histories still
  stop recovery. The reviewed files are hashed, backed up, and checked again;
  intermediate continuation records are imported and the owning chat is read
  through the isolated app-server. This is explicit acknowledgement of a reviewed
  chain, not a general force-import option. Retain the flag on subsequent pulls
  while a source still indexes an older segment of that chat.
- **New schema, changed lineage, project membership, missing ancestor, or failed
  app-server read:** keep the artifacts and investigate that explicit condition.
  The dated October 21 recovery scripts provide implementation history, not a
  reason to remove a failing check from the reusable utility.
- **Validation returned no turns:** inspect both the saved history rows and the
  rollout. A paginated index with no cached rows can return an empty history even
  when the transcript contains messages. Do not accept that as a successful
  recovery. An unrelated existing chat can be left untouched with
  `--exclude-thread THREAD_ID` while other chats are retrieved; report the
  exclusion and retain its transcript and database backups for separate repair.
- **Process interrupted during apply:** check `receipt.json`. An
  `apply-started.json` without a successful receipt means completion is unknown.
  Inspect both databases and the saved plan before continuing. Do not copy backup
  files over live SQLite/WAL files. A whole-store restoration requires stopping
  the affected backend, preserving its current databases, and reviewing changes
  since the backup; otherwise it could erase newer work.

SQLite transactions roll back handled errors, but a machine failure during the
commit of two databases in WAL mode is not a guaranteed all-or-nothing commit
across both files. The before-apply backup and completion receipt are therefore
part of the recovery procedure.

If a source's node-local database has already been erased, a startup pull cannot
recover that erased data. Existing durable snapshots can still be used explicitly:

```bash
python3.12 scripts/retrieve_codex_chats.py \
  --from-snapshot n0001=/absolute/path/to/saved/n0001 \
  --from-snapshot n0003=/absolute/path/to/saved/n0003
```

Each directory must contain both database files. Add `--apply` after inspecting
the result. Snapshot mode does not contact those nodes and makes no claim to
include messages written after the snapshots. It still needs the shared rollouts
and ancestors; if those are also missing, restore them separately from reviewed
backups before running this utility.

Pulling at startup removes the dependency on a clean session shutdown. Periodic
durable snapshots during long sessions would address the separate risk of losing
node-local storage before the next pull; no periodic job is installed here.

## Validation and maintenance

Run the data-preservation tests with:

```bash
python3.12 -m unittest discover -s tests/tools -v
```

They cover content conflicts, split histories, continuation-only rows, active
session exclusions, schema changes, concurrent destination changes and rollback,
organization preservation, repeat imports, WAL backups, writer locks, and reviewed
continuation chains (including missing ancestors, cycles, transcript changes, and
attempts to select an older or incomplete source).
Before changing merge behavior, also replay it on isolated copies of real saved
databases and rollouts and validate the resulting chats with the installed Codex
engine. Keep all live database paths out of that test.

Validated on October 21, 2026: all 12 tests passed. A full replay against isolated
copies of that day's snapshots imported 12 new records (11 user chats and one
helper thread), six history updates, and six metadata updates; all 24 chats
loaded through app-server. Repeating the command made no changes. An SSH preview
against n0001 and n0003 also passed without writing live chat records.
