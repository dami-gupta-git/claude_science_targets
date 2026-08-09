# Claude Science skills for drug target investigation

Custom Claude Science skills for judging whether a gene or compound is worth pursuing as a drug target.

Most analysis skills share one premise: the obvious analysis gives the wrong answer until a specific control is run. The skills encode those controls, the evidence order that puts cheap disconfirming data first, and the helpers needed to run them reproducibly inside Claude Science.

`skills/` holds the skill source (one folder per skill).  
`results/` holds analysis outputs, laid out `results/<topic>/<run>/`.

See [`skills/README_SKILLS.md`](skills/README_SKILLS.md) for a plain-language description of each skill, and [`skills/coding-standards`](skills/coding-standards) for the layout and naming rules both directories follow.

## Skills

Grouped as in `README_SKILLS.md`:

- **Judging a drug target** — `target-triage-public-data`, `depmap-local`, `opentargets-evidence`, `depmap-fusion`, `marker-contrast-null`, `depmap-genetea`, `target-lit-brief`
- **Molecules and structures** — `get-protein-structure`, `boltz-affinity-triage`
- **Evidence and claims** — `ai-origination-audit`
- **Conventions** — `coding-standards`, `doc-style`, `code-review`
- **Personal** — `morning`

Each skill directory contains:

- `SKILL.md` — instructions the agent loads
- `README.md` — human-readable documentation
- `kernel.py` (where present) — helpers preloaded into the Python kernel
- `tests/` (where present) — offline test suite

This repository is the authoritative source. `host.skills.publish()` uploads a directory into the skill registry; the runtime cache under `~/.claude-science/orgs/<org>/skills/` is a deployment target regenerated from it.

## Quick start on another Claude Science account

A git clone alone does not make the skills usable. Load them into the target account’s skill registry, then satisfy any data or connector requirements the skills you want actually need.

### 1. Publish each skill

In a `repl` cell, for each skill directory:

```python
for f in ["SKILL.md", "README.md", "kernel.py"]:  # kernel.py where present
    host.skills.edit("<skill-name>", f, open(f"skills/<skill-name>/{f}").read())
host.skills.publish("<skill-name>")
```

This is what makes `skill("<skill-name>")` resolve.

### 2. Install the Python environment (for tests and local runs)

```bash
uv sync                  # core test suites
uv sync --extra genetea  # + GeneTEA path for depmap-genetea
uv sync --all-extras     # + matplotlib for re-running scripts under results/
```

Verified on Python 3.12. On a bare checkout (no data files), 289 tests pass and 29 skip. With `$DEPMAP_ROOT` and `$GENETEA_MODEL` set and the `genetea` extra installed, `depmap-genetea` runs fully.

Working inside Claude Science rather than a clone, add packages with `manage_packages` instead of bare `pip install` (see `coding-standards`).

### 3. Point DepMap-dependent skills at a local release

`depmap-local` and `depmap-genetea` read CRISPR gene-effect and cell-line metadata from disk. DepMap cannot be downloaded programmatically (bot wall), so the release must be obtained manually.

Copy `env.example` to `.env` (or export the variables) and set:

- `$DEPMAP_ROOT` — directory containing `CRISPRGeneEffect.csv`, `Model.csv`, and the other files listed in each skill’s README
- `$GENETEA_MODEL` — path to the trained GeneTEA model (≈1.11 GB, conventionally `$DEPMAP_ROOT/genetea/GeneTEA.pkl`), from [Figshare 10.6084/m9.figshare.28635317](https://doi.org/10.6084/m9.figshare.28635317)

Both helpers also accept an explicit `root=` / `path=` argument. They raise `FileNotFoundError` naming the missing variable when unset; there is no machine-specific fallback.

### 4. Point results writers at this repo (optional but recommended)

Several skills write full runs (tables + README) under `results/<topic>/<run>/`. Set `$SCIENCE_RESULTS_ROOT` to this checkout’s `results/` directory, or pass `root=` on each call. See `env.example`.

### 5. Attach connectors

| Connector            | Needed by                |
|----------------------|--------------------------|
| `clinical-genomics`  | `opentargets-evidence`   |
| `chembl`             | `boltz-affinity-triage`  |
| `chemistry`          | `boltz-affinity-triage`  |

Authorize under Settings → Connectors. Skills reach them via `host.mcp(...)`. Nothing further to install.

`depmap-fusion` takes pre-computed statistics as arguments and needs no connector. `target-triage-public-data`’s Open Targets fallback hits the public GraphQL endpoint directly.

## What each skill needs

| Skill                        | Local DepMap | Connector              | Extra packages     | Live domains it calls                                      |
|------------------------------|--------------|------------------------|--------------------|------------------------------------------------------------|
| `depmap-local`               | required     | —                      | pyarrow            | —                                                          |
| `depmap-genetea`             | required + GeneTEA.pkl | —               | scikit-learn       | —                                                          |
| `depmap-fusion`              | —            | —                      | —                  | —                                                          |
| `target-triage-public-data`  | optional     | —                      | statsmodels        | api.platform.opentargets.org, ftp.sanger.ac.uk, cbioportal, GDC |
| `opentargets-evidence`       | —            | clinical-genomics      | —                  | (via connector)                                            |
| `marker-contrast-null`       | —            | —                      | —                  | —                                                          |
| `get-protein-structure`      | —            | —                      | gemmi              | RCSB, UniProt, EBI, AlphaFold DB                           |
| `boltz-affinity-triage`      | —            | chembl, chemistry      | rdkit              | (via connectors)                                           |
| `ai-origination-audit`       | —            | —                      | —                  | clinicaltrials.gov                                         |
| `target-lit-brief`           | —            | pubmed (control plane) | —                  | (via connector)                                            |
| `code-review`, `doc-style`, `coding-standards` | — | —                | —                  | —                                                          |
| `morning`                    | —            | calendar/email/chat    | —                  | (via connector)                                            |

Four skills need nothing beyond publishing: `code-review`, `doc-style`, `coding-standards`, and `morning` (the last needs a connector only to be useful). `marker-contrast-null` runs on the default environment with synthetic data.

## Tests

Seven skills ship test suites. Run from the skill directory:

```bash
python -m pytest
```

Each carries a `pytest.ini` that pins `rootdir`. Suites that need missing data files or models skip rather than fail, so a bare checkout stays green.

| Skill                       | Tests                  | Interpreter needs  |
|-----------------------------|------------------------|--------------------|
| `boltz-affinity-triage`     | 130 passed, 8 skipped  | rdkit              |
| `target-triage-public-data` | 24                     | scipy, statsmodels |
| `depmap-local`              | 19                     | pandas, pyarrow    |
| `depmap-fusion`             | 19                     | pandas             |
| `depmap-genetea`            | 13 passed, 27 skipped  | scikit-learn       |
| `marker-contrast-null`      | 16                     | scipy              |
| `get-protein-structure`     | 11                     | gemmi              |

## Scope and layout

- `skills/` contains only skill source. Do not leave analysis outputs beside a `SKILL.md`; they would ship to every machine that syncs the catalog.
- Analysis outputs belong under `results/<topic>/<run>/` (tables, figures, and the scripts that produced the run), per `coding-standards`.
- Publishing, agent profiles, and connector configuration are handled by the `customize` and `skill-creator` skills, not from this repository.

## Portability notes
Paths are driven by environment variables (`$DEPMAP_ROOT`, `$GENETEA_MODEL`, `$SCIENCE_RESULTS_ROOT`) or explicit arguments, with clear errors when unset. See `env.example` and `pyproject.toml`.
