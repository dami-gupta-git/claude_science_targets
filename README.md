# claude_science

Source for the custom Claude Science skills and analysis results built on this
account. `skills/` holds skill source, one folder per skill; `results/` holds
what each analysis produced, laid out `results/<topic>/<run>/`. See
`skills/README_SKILLS.md` for what each skill does, and `skills/coding-standards`
for the layout and naming rules both directories follow. This repository was
moved from `claude_science/` to `claude_science/claude_science_targets/`; any
doc below written before that move and still saying `claude_science/skills/…`
means this directory's `skills/…`.

## The `skills/` directory

Source for thirteen custom Claude Science skills, one directory per skill. Each
holds a `SKILL.md` that instructs the agent, a `README.md` for the person
reading the directory, and — where the skill ships helpers — a `kernel.py`
preloaded into the Python kernel plus a `tests/` suite. This is the authoritative
copy: `host.skills.publish()` uploads a directory to the skill registry, and the
runtime's cached copy under `~/.claude-science/orgs/<org>/skills/` is a
deployment target regenerated from it, not a source tree.

Most of the analysis skills answer a version of one question — is this gene or
compound worth pursuing — and share a premise: the obvious analysis gives the
wrong answer until a specific control is run. Three carry writing, coding and
review conventions rather than analysis, and one is a personal daily routine.

### Contents

- `README_SKILLS.md` — what each of the thirteen skills does, in plain language,
  grouped by purpose. Start here.
- `pyproject.toml`, `uv.lock`, `.python-version` — the environment the suites
  run under. Dependency comments name the skill that needs each package.
- `main.py` — `uv init` placeholder; no skill imports it.

Skills grouped as `README_SKILLS.md` groups them:

- Judging a drug target — `target-triage-public-data`, `depmap-local`,
  `opentargets-evidence`, `depmap-fusion`, `marker-contrast-null`,
  `depmap-genetea`.
- Molecules and structures — `get-protein-structure`, `boltz-affinity-triage`.
- Evidence and claims — `ai-origination-audit`.
- Conventions — `coding-standards`, `doc-style`, `code-review`.
- Personal — `morning`.

### Tests

Seven skills ship suites, run from the skill directory with `python -m pytest`.
Each carries a `pytest.ini` pinning `rootdir`; without one, pytest's upward
discovery reads this directory's `pyproject.toml` instead.

| Skill | Tests | Interpreter needs |
|---|---|---|
| `boltz-affinity-triage` | 130 passed, 8 skipped | rdkit |
| `target-triage-public-data` | 24 | scipy, statsmodels |
| `depmap-local` | 19 | pandas, pyarrow |
| `depmap-fusion` | 19 | pandas |
| `depmap-genetea` | 13 passed, 27 skipped | scikit-learn |
| `marker-contrast-null` | 16 | scipy |
| `get-protein-structure` | 11 | gemmi |

Suites needing a data file or a trained model skip rather than fail when it is
absent, so a bare checkout runs green. Setting `$DEPMAP_ROOT` and
`$GENETEA_MODEL` turns the 27 `depmap-genetea` skips into passes.

### Data locations

No skill carries a built-in data path. `depmap-local` and `depmap-genetea` read
a manually downloaded DepMap release via `$DEPMAP_ROOT`, and `depmap-genetea`
additionally needs a trained GeneTEA model via `$GENETEA_MODEL`; both raise
`FileNotFoundError` naming the variable when it is unset rather than trying a
path from the machine they were written on. See "Point `depmap-local` and
`depmap-genetea` at a DepMap release" below for `env.example` and the
per-skill requirement table.

### Scope

`skills/` holds skill source only. Analysis outputs — tables, figures and
the scripts wiring a particular run — belong in `results/<topic>/<run>/`, per
`coding-standards`; a file left beside a `SKILL.md` ships to every machine that
syncs the catalog. Publishing to the registry, agent profiles and connector
configuration are handled by the `customize` and `skill-creator` skills, not
from here.

## Running these on another Claude Science account

A checkout alone does not make a skill usable — `SKILL.md` in a cloned folder
is inert until it is loaded into the target account's skill registry. Four
things to do, then a per-skill requirement table for what each one needs beyond
that.

### 1. Publish each skill

In a `repl` cell, per skill directory:

```python
for f in ["SKILL.md", "README.md", "kernel.py"]:  # kernel.py where present
    host.skills.edit("<skill-name>", f, open(f"skills/<skill-name>/{f}").read())
host.skills.publish("<skill-name>")
```

This is what makes `skill("<skill-name>")` resolve; nothing else does.

### 2. Install the Python environment

`pyproject.toml` and `uv.lock` describe everything the test suites import, and
`.python-version` pins the interpreter the suites are verified on:

```bash
uv sync                  # all seven suites run; the GeneTEA tests skip
uv sync --extra genetea  # + the GeneTEA path, 27 further tests in depmap-genetea
uv sync --all-extras     # + matplotlib, for re-running the scripts under results/
```

The `genetea` extra installs GeneTEA from the Broad's repository at a pinned
commit, because it is not published on PyPI. The pin is load-bearing: the
trained model is a pickle written by that build. That extra also carries the five
libraries GeneTEA's `__init__` imports, since its own packaging declares no
dependencies and `import GeneTEA` fails without them.

Each dependency is there for one skill: `pyarrow` for `depmap-local`'s parquet
cache, `scikit-learn` for `depmap-genetea`, `rdkit` for `boltz-affinity-triage`,
`gemmi` for `get-protein-structure`, `statsmodels` for
`target-triage-public-data`. Working inside Claude Science rather than from a
clone, add them with `manage_packages` instead of a bare `pip install` — see
`coding-standards` for why.

Verified on Python 3.12: 289 tests pass and 29 skip on a bare checkout with no
data files present; with `$DEPMAP_ROOT` and `$GENETEA_MODEL` set and the
`genetea` extra installed, `depmap-genetea` runs 41 passed and 0 skipped.

### 3. Point `depmap-local` and `depmap-genetea` at a DepMap release

Both read CRISPR gene effect and cell-line metadata from files on disk, not
over the network — `depmap.org` sits behind a bot-verification wall, so the
release has to be downloaded manually and cannot be scripted. Copy
`env.example` to `.env` (or export directly) and set `$DEPMAP_ROOT` to the
directory holding `CRISPRGeneEffect.csv`, `Model.csv` and the other files each
skill's README lists; both skills also accept a `root=` argument as an
override. Neither has a built-in default — `depmap_root()` raises
`FileNotFoundError` naming the variable when it is unset, rather than trying a
path from the machine these were developed on.

`depmap-genetea` additionally needs the trained GeneTEA model — 1.11 GB, at
`$GENETEA_MODEL` (conventionally `$DEPMAP_ROOT/genetea/GeneTEA.pkl`), from
Figshare
[10.6084/m9.figshare.28635317](https://doi.org/10.6084/m9.figshare.28635317).
Same behaviour: `genetea_load()` raises naming `$GENETEA_MODEL` when neither it
nor `path=` is supplied.

### 4. Point every skill's results-writer at this repo

Seven skills (`target-lit-brief`, `boltz-affinity-triage`, `depmap-fusion`,
`depmap-genetea`, `get-protein-structure`, `marker-contrast-null`,
`ai-origination-audit`) save a full run — tables plus a README — into this
repo's `results/<topic>/<run>/` via a `*_run_dir()` helper. Set
`$SCIENCE_RESULTS_ROOT` to this checkout's `results/` directory (see
`env.example`), or pass `root=` explicitly on each call. No cwd-relative
default: each raises `FileNotFoundError` naming the variable when neither is
set, rather than silently creating a new `results/` folder wherever the
kernel session's working directory happens to be.

### 5. Attach connectors

Settings → Connectors, for whichever skills you intend to use:

| Connector | Needed by |
|---|---|
| `clinical-genomics` | `opentargets-evidence` |
| `chembl` | `boltz-affinity-triage` |
| `chemistry` | `boltz-affinity-triage` |

Skills reach these via `host.mcp(...)` in the `repl` tool; there is nothing to
install for a connector, only to authorize.

`depmap-fusion` calls no connector — `fuse_target_row()` takes evidence and
dependency statistics as arguments, from whichever source produced them.
`target-triage-public-data`'s Open Targets fallback (see step 2 above) hits
`api.platform.opentargets.org` directly with `urllib.request.urlopen`, not
through a connector.

## What each skill needs

Derived from each skill's own imports and network calls, not from its prose —
grep against `kernel.py` undercounts when a requirement is documented but not
called from code, and overcounts when a domain name appears only in an example.

| Skill | Local DepMap release | Connector | Extra packages | Live domains it calls |
|---|---|---|---|---|
| `depmap-local` | **required** | — | `pyarrow` | — |
| `depmap-genetea` | **required**, plus `GeneTEA.pkl` | — | `scikit-learn` | — |
| `depmap-fusion` | — (takes stats as a dict argument; source-agnostic) | — | — | — |
| `target-triage-public-data` | optional — prefers it, falls back to Open Targets GraphQL for the same scores at reduced line coverage | — (fallback hits `api.platform.opentargets.org` directly via `urllib`, not a connector) | `statsmodels` | `api.platform.opentargets.org`, `ftp.sanger.ac.uk` (GDSC), `www.cbioportal.org`, `api.gdc.cancer.gov` |
| `opentargets-evidence` | — | `clinical-genomics` | — | (via connector) |
| `marker-contrast-null` | — | — | — | — |
| `get-protein-structure` | — | — | `gemmi` | `data.rcsb.org`, `files.rcsb.org`, `rest.uniprot.org`, `www.ebi.ac.uk`, `alphafold.ebi.ac.uk` |
| `boltz-affinity-triage` | — | `chembl`, `chemistry` | `rdkit` | (via connectors) |
| `ai-origination-audit` | — | — | — | `clinicaltrials.gov` |
| `code-review` | — | — | — | — |
| `doc-style` | — | — | — | — |
| `coding-standards` | — | — | — | — |
| `morning` | — | calendar/email/chat connector (whichever is attached) | — | (via connector) |

Four skills run with nothing beyond publishing: `code-review`, `doc-style`,
`coding-standards`, `morning` (the last needs a connector only to do anything
useful, not to load). `marker-contrast-null` needs only `numpy`/`pandas`/`scipy`,
already in the default environment — its own test suite runs on synthetic data
with no external dependency at all.

## Portability fixes applied

Three absolute-path defaults tied this repository to the machine it was
developed on; all three are now environment-driven with no machine-specific
fallback:

- `depmap-local/kernel.py`'s `depmap_root()` and `depmap-genetea/kernel.py`'s
  own copy of the same function (a separate implementation, not a shared
  import) both required `$DEPMAP_ROOT` or `root=`, raising `FileNotFoundError`
  naming the variable otherwise.
- `depmap-genetea/kernel.py`'s `genetea_load()` requires `$GENETEA_MODEL` or
  `path=` the same way.
- `coding-standards/SKILL.md`'s layout rule refers to "this repository's
  `skills/…`" rather than one developer's absolute path.

`env.example` and `pyproject.toml` at the repo root capture the variables and
packages these fixes now require explicitly.
