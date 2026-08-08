---
name: doc-style
description: House writing conventions for every written deliverable — README, SKILL.md, report, results summary, dossier, methods section, docstring, changelog, or documentation of any kind. Specifies document structure (brief plain-prose overview paragraph, then bullet points naming specific functions or files), professional non-opinionated tone, and what to leave out (bug changelogs, arrow-and-fragment shorthand, marketing adjectives, self-narration). Load BEFORE writing or revising any prose deliverable, README, or documentation file, and before summarising results in a saved report — including when the user says write a readme, document this, write it up, summarise the results, or add docs.
---

# Documentation style

House conventions for prose deliverables. Apply to READMEs, skill
documentation, result summaries, and reports.

## Structure

Open with a **short overview paragraph** — three to five sentences of plain
prose stating what the thing does and the problem it addresses. Not a feature
list, not a chain of semicolons, and not a single fragment.

Follow with **bullet points naming specific functions or files**, each with a
one-line description of what it does:

```
- `get_structure(query)` — identifier to prepared file: resolves the protein,
  ranks all mapped entries, downloads and cleans the best.
- `rank_structures(accession)` — the ranked candidate table alone.
```

Name the actual callable or filename, not the capability in the abstract.

Keep the whole document proportionate to what it documents. A README for a
single module should be readable in a couple of minutes.

## Tone

Professional and descriptive. State what the code does and what the data
shows; leave the reader to judge significance.

Avoid:

- Editorialising and value judgements — "the interesting part", "worth knowing
  about", "the row that settles it", "a tempting substitute".
- Self-narration about the development process — "an earlier pass had this
  wrong", "found in review", "I rejected this approach".
- Marketing register — "powerful", "seamless", "robust" as unquantified
  adjectives.
- Emphasis by adjective where a number would do. "13.5 A RMSD and 0.18
  recovery" carries more than "poor agreement".

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

Where a threshold or grade is reported, say what calibrated it and how to
re-derive it.

## Scope statement

Close with what the thing does **not** do, naming the tools that do handle it.
One short paragraph, factual, no apology.
