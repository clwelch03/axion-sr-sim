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
   alone do not trigger imports or history conflicts. Updated/recency timestamp
   differences alone also leave an existing chat untouched; those timestamps are
   imported when its content or other metadata changes. Review/helper threads
   are included in counts.
   If a shared rollout has already migrated to paginated history but a node's
   index still says `legacy`, the plan records a one-way index repair. The
   rollout header must identify that chat and declare `paginated`; its hash is
   backed up and checked again before applying. Genuine legacy rollouts remain
   legacy. This also repairs an affected local index without rewriting rollouts.
   New project/section memberships require manual review. Different continuation
   paths are resolved automatically when one eligible source descends from every
   known rollout of the same chat through `history_base`. The script checks the
   ancestor identities and byte/ordinal boundaries, rather than choosing a path
   by timestamp. Unrelated or competing branches still stop the import.
   A continuation can branch before the end of an older segment. The terminal
   and `plan.json` report those older records and turn IDs; the entire ancestor
   rollout is backed up and all known database rows are retained. The chat opens
   with the selected source's continuation history, so an older tail outside its
   boundary may not appear in that view. The script does not splice it into the
   continuation or rewrite shared rollouts.
4. Backs up affected rollouts and their ancestors with SHA-256 hashes. It merges
   into temporary database copies and verifies the exact imported rows. An
   isolated Codex app-server reads every affected chat, including the owners of
   continuation-only history. If an imported paginated chat has no cached history
   or projection state, and its rollout is a complete original segment, the
   engine rebuilds the cache by resuming a disposable copy. Only rows within the
   original transcript's byte/ordinal bounds are retained. Any settings events
   appended to the copy are discarded, and the rebuilt cache is validated again
   against a fresh copy of the original transcript. Partial caches and
   continuation segments are not rebuilt this way. Legacy chats also require
   resuming a copy before this engine can read their turns.
   Validation and cache rebuilding do not start model turns or resume live chats.
   Each probe runs from its disposable directory in a separate process session.
   Shutdown also stops its background helpers, including plugin-clone processes,
   before scratch cleanup. A scratch directory that cannot be removed is reported
   with a warning and its retained path; it does not invalidate successful history
   checks or hide a validation failure.
   A zero-turn
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
The installed engine supports `thread/read` with `includeTurns: true` for
paginated chats; legacy chats need an isolated `thread/resume` first. Current
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
| `plan.json` | Proposed changes, source labels, excluded IDs, conflicts, verified continuation hashes and boundary reports |
| `transcripts/`, `transcripts.json` | Affected rollouts, ancestors, and hashes |
| `validation.json`, `validation.stderr.log` | Isolated app-server results and diagnostics |
| `rebuild-validation.json`, `rebuild-validation.stderr.log` | Isolated cache rebuild results, when needed; recovered rows are recorded in `plan.json` |
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
- **Different rollout/continuation paths:** rerun the normal command with the
  corrected script. A stale node indexing an earlier segment of a verified chain
  no longer requires a flag, including on repeated pulls. If verification fails,
  the terminal prints the conflicting title, ID, and reason; `plan.json` also
  lists each source's path. Review the shared rollout headers and their
  `history_base` links, including the content around continuation boundaries.
  For an explicitly reviewed chain whose boundary format the automatic check
  cannot verify, select that specific chat and source with
  `--reviewed-continuation THREAD_ID=SOURCE` (repeatable), using the source label
  from `plan.json`, such as `n0003_scs00`. Add this to the normal command with
  `--apply`. The utility still requires the selected source to contain all known
  content rows and verifies that its chain includes every older rollout for the
  same conversation. Missing, cyclic, unrelated, or divergent histories still
  stop recovery. The reviewed files are hashed, backed up, and checked again;
  intermediate continuation records are imported and the owning chat is read
  through the isolated app-server. This is explicit acknowledgement of a reviewed
  chain, not a general force-import option. It bypasses automatic byte/ordinal
  boundary checks, but still verifies same-chat ancestry and content preservation.
- **New schema, changed lineage, project membership, missing ancestor, or failed
  app-server read:** keep the artifacts and investigate that explicit condition.
  The dated October 21 recovery scripts provide implementation history, not a
  reason to remove a failing check from the reusable utility.
- **Validation returned no turns:** inspect both the saved history rows and the
  rollout. A paginated index with no cached rows can return an empty history even
  when the transcript contains messages. The script now rebuilds a wholly absent
  cache for an imported original segment using the installed engine on isolated
  copies, then validates it again. Other empty histories still stop recovery;
  do not accept those as a successful import. An unrelated existing chat can be left untouched with
  `--exclude-thread THREAD_ID` while other chats are retrieved; report the
  exclusion and retain its transcript and database backups for separate repair.
- **`thread/resume: ... list_turns is not supported yet`:** a node can still
  index a chat as `legacy` after another node migrates its shared rollout to
  paginated history. Use the corrected script and rerun the same command. It
  verifies the rollout header, repairs the index flag, and validates the full
  history on disposable copies before importing. The failed validation occurs
  before live imports; do not delete chat records or migrate shared rollouts
  to work around it.
- **`Directory not empty` under `/tmp/codex-chat-probe-*/home/.tmp/plugins-clone-*`:**
  older versions waited for app-server alone, leaving background plugin-clone
  processes writing during temporary-directory deletion. Use the corrected
  script and rerun the normal command. This failure happens after validation but
  before live import. Probe shutdown now stops its whole private process group,
  retries scratch cleanup, and warns if disposable files must be retained.
  Validation, source rechecks, writer locks, backups, and import verification still
  run; a scratch warning does not bypass any of them.
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
organization preservation, repeat imports, WAL backups, writer locks, and automatic
and reviewed continuation chains (including parent tails, invalid boundaries,
missing ancestors, cycles, transcript changes, and attempts to select an older
or incomplete source).
They also cover timestamp-only differences, missing-cache recovery, preservation
of source storage, and rejection of rebuilt rows outside the original transcript.
Before changing merge behavior, also replay it on isolated copies of real saved
databases and rollouts and validate the resulting chats with the installed Codex
engine. Keep all live database paths out of that test.

Validated on October 21, 2026: all 12 tests passed. A full replay against isolated
copies of that day's snapshots imported 12 new records (11 user chats and one
helper thread), six history updates, and six metadata updates; all 24 chats
loaded through app-server. Repeating the command made no changes. An SSH preview
against n0001 and n0003 also passed without writing live chat records.

Validated on October 29, 2026 against isolated copies of the failed
`20261029T163304Z-_fkw5c5w` snapshots: imported 13 new chats, eight updates, and
two continuation records without a reviewed-continuation flag. All 21 affected
chats loaded through the installed engine, including the rebuilt guardian review.
The older turn outside the radial BVP continuation boundary remained in the
destination history tables. A second pull made no changes. Live stores and
shared rollouts were not modified by this replay.

Validated on November 2, 2026: all 55 recovery tests passed, including real
subprocess checks for background writers after successful and failed validation,
and a complete import/repeat with simulated SSH and app-server responses while
scratch deletion was forced to fail. Replaying the failed
`20261102T172824Z-f2b9b0te` snapshots through the installed engine returned the
same 156 chats and 1,323 turns, with identical turn IDs and successful scratch
cleanup. The normal live `--apply` command on n0003 then completed with 40 new
chats, 116 updates, and eight continuation records, after fresh snapshots and
rechecks of n0001 and n0002. Its backups and successful receipt are in
`recovery/chat-retrieval/20261102T173455Z-gtfte2ue/`. Repeating the same command
reported `Already up to date; no live databases changed`, with the repeat receipt
in `recovery/chat-retrieval/20261102T173530Z-0aahsuif/`.
