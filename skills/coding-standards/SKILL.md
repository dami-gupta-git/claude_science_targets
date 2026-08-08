---
name: coding-standards
description: House coding conventions for analysis code in Claude Science — where source files must be written, how to structure a reusable script, and what has to be true before an analysis is called finished. The central rule is that code lives on disk, not only in a kernel: any script, module, test, or helper worth running twice is written to a file and saved, never left as cell-only text that dies with the kernel. Also requires that every skill ship a README.md alongside its SKILL.md, written to the doc-style conventions. Load BEFORE writing analysis code, building a pipeline, authoring or publishing a skill, writing a module or test suite, or refactoring a notebook-style cell into something reusable — and whenever the user says write a script, make this reproducible, save the code, turn this into a module, add tests, write a skill, or asks why an earlier analysis cannot be re-run.
---

# Coding standards

Analysis code in this environment is written twice by default: once into a
kernel cell so it runs now, and once onto disk so it can be run again. Only the
second one survives. A `python` cell's namespace is gone at kernel restart, and
the task workspace is swept after idle gaps — so a function that exists only in
a cell is a result you cannot reproduce and a colleague cannot inspect. Treat
"it ran and I saw the number" as the beginning of the work, not the end.

The failure is easy to spot after the fact. A project directory holding 31 CSVs,
7 figures and one `.py` file is a project whose outputs were all persisted and
whose code was almost entirely thrown away.

## Where code goes

Two destinations, and code that matters goes to both:

1. **The artifact store**, via `save_artifacts(files=[...], language=...)`, so
   the file is versioned, visible to the user, and linked to the session that
   produced it. Update in place with
   `version_of={"script.py": "<artifact_id>"}` instead of accumulating
   `script_v2.py`, `script_final.py`.
2. **A granted host path**, for a copy in the user's real filesystem — openable
   in an editor, committable, usable outside this app. Check
   `list_host_grants()` first, write with `edit_file` (atomic replace, original
   preserved on failure), and follow the layout already in the directory rather
   than inventing one.

Order matters, because `save_artifacts` resolves paths against the workspace or
a registered repo root only: an absolute path into a granted host directory is
rejected with *"not under an allowed root"*. So write the file in the workspace,
save it by its relative name, then copy it to the host path — or register the
repo with `manage_environments(mode="register", ...)` if you want absolute paths
to resolve. Saving from a host path directly does not work.

Neither destination is the workspace alone. Writing `helper.py` with a relative
path and stopping there persists nothing.

Never `rm` inside a granted host path — use `delete_host_files`, which prompts
the user and moves to Trash. An rw grant lets `rm` succeed silently on the
user's machine with no undo.

## Skill source lives in the user's repo, not the registry

`host.skills.publish()` writes to a server-side registry and the runtime then
materialises a copy under `~/.claude-science/orgs/<org>/skills/<name>/`. That
cached copy carries `.sync-org` and `.catalog_stamp` markers, is regenerated
from the registry, and is not a source tree — its file timestamps reflect the
last catalog sync, not when the skill was written, so they cannot be used to
date authorship. A skill that exists only in those two places has no diffable,
reviewable, user-controlled copy anywhere.

So a skill is authored **into the user's own directory first**: this
repository's `skills/<skill-name>/` — the checkout holding this file, at
whatever path it lives on the current machine (e.g.
`~/code/portfolio/claude_science/skills/<skill-name>/`; confirm the actual path
from the project's granted host paths or by asking rather than assuming that
example). Write `SKILL.md`, `README.md`, `kernel.py` and `tests/` there, run the
suite from that directory, and publish from those files. The repo copy is the
source; the registry is the deployment target.

Exclude the runtime's own bookkeeping when copying anything out of the cache —
`.sync-org`, `.catalog_stamp`, `__pycache__`, `.pytest_cache` are per-machine
and only add diff noise.

**Do not run git.** No `init`, `add`, `commit`, `branch`, or `push` unless the
user asks in that turn. Write the files, say what changed, and leave staging and
history to them.

## Results go to `results/`, never inside a skill

This repository's root holds exactly two directories: `skills/` for skill
source and `results/` for everything an analysis produced.
Results are laid out `results/<topic>/<run>/` — topic names the analysis type,
run names the subject:

```
results/target_triage/usp1/          tables, figures and README for one run
results/target_triage/usp1/scripts/  the code that produced them
results/protein_structure/wrn/
```

Every run directory carries a `README.md`, and every topic directory carries a
short index naming each run and what it found. The run README opens with a
`Result` section stating the finding in prose — that section is the deliverable;
the tables are its evidence — then `Files` with a line per file, `Data sources`
naming the skills and data release used, and `Limits` for what the result does
not cover. Use `README.md`, not `RESULTS.md`: it is what renders when the folder
is opened, and `doc-style` already governs it.

Code that produced a run goes in that run's `scripts/` subfolder, so the run
directory itself reads as data. This is analysis-specific wiring — the gene
list, the marker panel, the indications tested — and it stays with its outputs
rather than moving into a skill. When part of it turns out to be reusable,
extract that part into a skill and leave the wiring behind: the USP1 run's two
controls became `marker-contrast-null` this way, while its comparator lists
stayed in `results/target_triage/usp1/scripts/`.

A second run of the same analysis becomes a sibling under the same topic
(`target_triage/dctpp1/`), so runs stay comparable and the topic folder does not
accumulate loose files.

Naming follows the two directory trees. Skills are lowercase-with-hyphens,
because the Agent Skills standard requires the directory name to match the
`name:` field in `SKILL.md` frontmatter. Result directories are `snake_case`,
which keeps them importable as Python packages and makes the separator alone say
which tree a path belongs to.

**Never write results into a skill directory.** `host.skills.publish()` uploads
the whole directory, so a CSV or figure left beside `SKILL.md` ships to the
registry and materialises on every machine that syncs the catalog, permanently.
The cardinality is also wrong in both directions: one skill produces many runs,
and one run uses many skills, so no single skill directory is the right home for
a given output. The exception is a small fixture that a documented threshold was
calibrated from — that belongs with the skill it calibrates, because both
`coding-standards` and `doc-style` require a threshold to state what set it.

Per-run provenance — a copy of the skill as it stood when the analysis ran —
does not go under `results/` either. Record the version in the run's README
instead; a second copy of a skill on disk is indistinguishable from a live one
and will eventually be edited or published by mistake.

## What earns a file

Not every line needs a module. The test is whether the code will be *run again*
— by a later cell, a later session, or the user:

- **Write a file** for anything reusable: a function called more than once, a
  parser, a plotting routine, a pipeline stage, a fetch-and-clean step, a test.
- **Leave it in the cell** for genuine one-offs: `df.describe()`, a quick
  distribution check, printing a shape to decide the next line.

When a cell-only helper turns out to be needed a second time, that is the signal
to move it to a file — not to paste it again. Copy-paste across cells is the
common way a project ends up with three subtly different versions of the same
function and no way to tell which produced the saved figure.

## Structure

Keep imports at the top of the module, group a file around one concern, and give
functions arguments rather than reaching for globals — a function that reads `df`
from the enclosing namespace cannot be called from anywhere else, which defeats
the point of the file.

Docstrings carry the *why*. A signature already says what the arguments are; the
docstring should say what the function is for and name the constraint that
shaped it — the qualifier that must not be dropped, the field that is absent on
most rows, the cap the upstream tool enforces. That is the knowledge a future
caller cannot reconstruct from the code.

Prefer explicit failure to a silent fallback. When a caller names something that
does not exist — a chain, a column, a sample index — return or raise an error
that names what was asked for and what is available, rather than substituting
the first available item and reporting success. A quietly wrong answer costs far
more than a loud one.

Validate inline as you go: `assert len(df) > 0, f"got {df.shape}"` costs nothing
and turns a silent empty-frame bug into an immediate stop.

## Correctness patterns that keep recurring

Kernel reviews across this repository keep finding the same three bug shapes in
unrelated skills, each already fixed once elsewhere in the same codebase. Check
for these explicitly when writing a function and when fixing an edge case in
one — not just at first-write time.

A `dict.get(key, default)` default only fires when the key is absent, not when
its value is `None`. `effect_stats.get("frac_dependent", 0.0)` looks like it
defaults missing data to `0.0`, but a row carrying `"frac_dependent": None`
sails through untouched, and every `>=`/`<` comparison on it downstream raises
`TypeError`. A field that can be `None` needs its own explicit coalesce
(`x = d.get(k); x = default if x is None else x`), not a `.get(k, default)`
that only guards the absent-key case.

Every comparator side of a statistical test needs the same size guard, not just
whichever side failed first. `scipy.stats.mannwhitneyu` on an empty or
near-empty sample returns `NaN` rather than raising. If that `NaN` reaches a
"floor tiny p-values so they stay reportable" helper
(`p if p > 0 else sys.float_info.min`), the comparison silently becomes the
*most significant possible result* instead of an error — `NaN > 0` is `False`
in Python, so the floor branch fires. When one branch of a two-sample
comparison is guarded (`if n_pos < min_n: return "untestable"`), check that
every other branch of the same comparison has the identical guard before
treating the function as done.

Boundary and split artifacts get proven on the happy path, not the seam. A
regex matched as a bare substring (`"generative"` inside `"regenerative"`)
instead of word-bounded; a sentence's trailing period consumed by a `"..."`
join separator and never reattached to the fragment that owned it. Both were
caught only by hand-tracing the actual split or match output on adjacent,
non-duplicate input — a test that normalises away punctuation or word
boundaries before comparing (`.strip(".; ")`, case-insensitive substring) can
pass while the underlying bug still ships.

After fixing one instance, grep the file for the same call or idiom before
moving on. In three of five reviewed skills, the correct guard already existed
a few functions away in the same file, written for the identical library call.
One function guarded its `mannwhitneyu` comparator against empty samples with
a comment explaining why; a sibling function a few functions later ran the
same comparison unguarded. The fix pattern was already in the file — it just
had not been checked against every other caller of the same function.

## A skill ships with a README

Authoring a skill means writing three things, not one: `SKILL.md` for the agent
that will load it, `kernel.py` for any helpers it needs, and `README.md` for the
human reading the directory. The README is not a duplicate of `SKILL.md` — they
have different readers and different jobs. `SKILL.md` instructs an agent
mid-task; `README.md` tells a person what this directory is, what each file in
it does, and what calibrated any threshold it reports.

Write it to the conventions in `doc-style`, which govern it as they govern any
prose deliverable: load that skill first and follow it. In short — open with a
short plain-prose overview paragraph rather than a feature list, then bullet
points naming the actual functions and files with one line each, keep the tone
descriptive and free of marketing adjectives and development narration, ground
every quantitative claim in a saved artifact, and close with a scope statement
naming the tools that handle what this one does not.

Two rules specific to a skill README. State what calibrated any threshold the
skill applies and how to re-derive it, since a number with no provenance cannot
be revisited when the data changes. And where the skill ships tests, name the
command that runs them and the interpreter they need — a suite whose environment
is undocumented stops being run.

Skills whose behaviour depends on helper functions should list them by name and
signature, so the README stays checkable: every helper named in it must exist in
`kernel.py`, and every helper in `kernel.py` should be named.

## Reproducibility

The kernel tracks lineage through namespace variables, which constrains how
plots and fetches are written:

- Save figures from the figure object — `fig.savefig(...)`, not
  `plt.savefig(...)`; in R, `ggsave("out.png", plot = p)`. Bare `plt.*` calls
  produce broken lineage.
- Put a network fetch in its own cell and read the file in the next, so the
  fetch can be stubbed on replay.
- Pin the environment: pass `environment=` on every call, and
  `environment=` to `save_artifacts` to capture the package snapshot alongside
  the file.

State setup once. Within an environment, variables and imports persist across
calls, so re-emitting `import pandas as pd` in every cell is noise; the next
cell should be the incremental step.

## Tests

Anything with a correctness condition worth stating gets a test file, saved like
any other source. A test suite that runs without a GPU, a network connection, or
credentials will still run months later, so prefer synthesised inputs shaped
like the real thing over fixtures that need a live service.

Passing tests are weak evidence on their own. Confirm a test constrains what it
claims by breaking the behaviour deliberately and checking that the expected
test — and only that test — fails. A test that still passes with the logic
inverted is documentation, not a check.

Test the boundaries and the degraded paths, not just the happy one: empty
inputs, a missing file, a value exactly at a threshold, an absent optional
dependency. Floating-point thresholds deserve explicit attention, since a
difference that should exactly meet a limit can land just below it.

## Before calling it finished

- The code that produced every reported number exists in a file, and that file
  is saved.
- Re-running that file from a clean kernel reproduces the result.
- Numbers quoted in prose were read back from the saved output, not retyped
  from memory.
- Paths in the file are relative or derived, not absolute paths into a
  workspace that will not exist next time.
- No credential, token, or key appears in the source.
- A skill being published has a `README.md` beside its `SKILL.md`, and every
  helper the README names exists in `kernel.py`.

## Scope

This skill covers where code lives and how it is structured; it is not a
language style guide and does not specify formatting, naming, or line length —
follow the conventions already present in the file being edited. This skill
requires a README with every skill but does not restate how to write one: the
prose conventions for READMEs, reports and any other written deliverable are
`doc-style`'s, and that skill should be loaded before writing one. Figure
correctness is governed by `figure-style`; skill packaging, publishing, and the
`kernel.py` sidecar restrictions by `skill-creator` and `customize`.
