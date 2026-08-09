---
name: doc-style
description: House writing conventions for every written deliverable — README, SKILL.md, report, results summary, dossier, methods section, docstring, changelog, or documentation of any kind. Specifies document structure (brief plain-prose overview paragraph, then a table naming specific functions or files), professional non-opinionated tone, and what to leave out (bug changelogs, arrow-and-fragment shorthand, marketing adjectives, self-narration). Load BEFORE writing or revising any prose deliverable, README, or documentation file, and before summarising results in a saved report — including when the user says write a readme, document this, write it up, summarise the results, or add docs.
---

# Documentation style

House conventions for prose deliverables. Apply to READMEs, skill
documentation, result summaries, and reports.

## Structure

Open with a **short overview paragraph** — three to five sentences of plain
prose stating what the thing does and the problem it addresses. Not a feature
list, not a chain of semicolons, and not a single fragment.

Follow with a **table naming specific functions or files**, one row each: the
callable or filename in one column, a one-line description of what it does in
the other.

README files should be crisp. NO nitty-gritty details.

```
| Function | Description |
| --- | --- |
| `get_structure(query)` | identifier to prepared file: resolves the protein, ranks all mapped entries, downloads and cleans the best. |
| `rank_structures(accession)` | the ranked candidate table alone. |
```

A table over a bullet list here because the pattern is the same two-field
record repeated — name, then what it does — and a table lets the reader scan
the names down one column without re-reading the prose. Reserve a bullet list
for entries whose description does not fit a single cell (multi-sentence
caveats, worked examples) rather than forcing it into a table row.

Name the actual callable or filename, not the capability in the abstract.

Keep the whole document proportionate to what it documents. A README for a
single module should be readable in a couple of minutes.

## IMPORTANT: Tone

Professional and descriptive. State what the code does and what the data
shows; leave the reader to judge significance.

Strict Rule: Avoid these:

- Editorialising and value judgements — "the interesting part", "worth knowing
  about", "the row that settles it", "a tempting substitute".
- Self-narration about the development process — "an earlier pass had this
  wrong", "found in review", "I rejected this approach".
- Marketing register — "powerful", "seamless", "robust" as unquantified
  adjectives.
- Emphasis by adjective where a number would do. "13.5 A RMSD and 0.18
  recovery" carries more than "poor agreement".
- Technical-sounding labels standing in for the property that makes them true —
  "the helpers here are pure", "a robust parser", "an idempotent writer". Name
  the actual property instead: "none of these call the network", "re-running
  it on the same input overwrites rather than appends". A label can be true
  and still read as a claim about the code's quality rather than a fact about
  its behaviour; the reader should not have to take the label on faith.

Prefer the passive or impersonal voice for design rationale: "a LIGSITE-style
detector was evaluated and not adopted" over "I tried it and it did not work".


## Omit

- **Bug changelogs.** A README documents current behaviour, not repair history.
  Where a past defect left a constraint worth preserving, state the constraint
  as a property of the data or library — "residue numbers are only unique
  within a chain" — not as a fixed bug.
- **Arrow-and-fragment shorthand** in prose. `Identifier -> usable file.
  Resolves X; ranks Y; writes Z.` is a feature list. Arrows are acceptable
  inside tables and code.
- **Restating the obvious layout.** A table listing `README.md | This file`
  earns its place only when the directory has non-obvious contents.

## Numbers and claims

Every quantitative claim in a document must be reproducible from a saved
artifact. Before publishing, cross-check each figure against the data file that
produced it; where a table and a figure disagree because of filtering, state
the reason rather than silently letting them differ.

Avoid statements like '5 helper functions', or '2 files'. These change frequently, do not
hardcode.

Where a threshold or grade is reported, say what calibrated it and how to
re-derive it.

Do not state a test count — "63 tests", "130 passed, 8 skipped" — in a README.
The suite changes size independently of the prose describing it, so the number
goes stale the next time a test is added or removed and nothing forces the
README to update with it. Point at the command that runs the suite instead
("run `python -m pytest` from the skill directory") and let the test output be
the source of truth for how many there are.

The same restraint applies to describing the suite's contents, not just its
size: a README does not need a per-test or per-mutation breakdown of what the
suite checks — which named test catches which regression, or a table of
inversions tried against a table of tests that failed. That level of detail
belongs in the test file itself, next to the assertions it describes, where it
cannot drift out of sync with a rename or a new test. State in one line what
kind of thing the suite covers and how to run it; let the test names and
docstrings carry the rest.

## Other rules
- When asked to make a document 'briefer' or 'shorten', do not convert sentences into phrases.
That just makes it illegible. Cut detail.

## Scope statement

Close with what the thing does **not** do, naming the tools that do handle it.
One short paragraph, factual, no apology.
