# coding-standards

House conventions for analysis code written in Claude Science. The environment
makes it easy to produce a correct result that cannot be reproduced: a kernel
cell runs, prints a number, and the function that computed it disappears when the
kernel restarts or the task workspace is swept. This skill states where code must
be written so that does not happen, what structure a reusable file should have,
and what has to hold before an analysis is treated as finished.

The rule the rest follows from is that code lives on disk rather than only in a
kernel. Source files are written to a granted host path so they exist in the
user's own filesystem, and saved to the artifact store so they are versioned and
attributable to the session that produced them.

## Contents

- `SKILL.md` — the conventions: destinations for source files, which code earns a
  file, module structure, the reproducibility constraints imposed by lineage
  tracking, testing expectations, and a pre-completion checklist.

## Scope of each section

- **Where code goes** — `save_artifacts` by relative workspace path for the
  versioned copy, then `list_host_grants()` and `edit_file` for the host copy,
  `version_of=` for updates in place, and `delete_host_files` rather than `rm`
  under a grant. The order is load-bearing: `save_artifacts` resolves paths
  against the workspace or a registered repo root, so an absolute path into a
  granted host directory is refused.
- **A skill ships with a README** — `SKILL.md`, `kernel.py` and `README.md` have
  different readers; the README additionally records what calibrated each
  threshold and how to run the tests. Prose conventions are delegated to
  `doc-style` rather than restated.
- **What earns a file** — the reuse test that separates a genuine one-off cell
  from code that should be a module.
- **Structure** — argument-passing over namespace capture, docstrings that record
  the constraint rather than restate the signature, explicit errors over silent
  fallbacks, and inline assertions.
- **Reproducibility** — `fig.savefig` over `plt.savefig`, isolated fetch cells,
  and environment pinning.
- **Tests** — offline-runnable suites, boundary and degraded-path coverage, and
  confirming a test constrains behaviour by breaking that behaviour deliberately.
- **Before calling it finished** — the conditions to check before reporting a
  result.

## Basis for the README rule

Of the 14 editable skills installed when this convention was written, 7 shipped a
`README.md` and 7 did not. The seven with one were the skills carrying helper
modules and test suites; the omission fell on skills where a reader most needs to
know which file does what. Requiring the README at authoring time removes the
question of whether a given skill happens to have one.

## Basis for the persistence rule

The convention is not stylistic. In the project directory this skill was written
against, a file census across five completed analyses found 31 CSVs, 10 markdown
files and 7 figures, against a single `.py` file. Outputs were persisted in every
case; the code that generated them was retained once. Each of those analyses is
therefore readable but not re-runnable, and a reported number cannot be traced to
the operation that produced it.

## Scope

This skill does not specify formatting, naming, or line length, and is not a
language style guide — existing conventions in the file being edited take
precedence. Prose deliverables are governed by `doc-style`, figure correctness by
`figure-style`, and skill packaging including `kernel.py` sidecar restrictions by
`skill-creator` and `customize`.
