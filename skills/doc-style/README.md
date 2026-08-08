# doc-style

House writing conventions for prose deliverables in this repository — READMEs,
`SKILL.md` files, run reports, result summaries and docstrings. The conventions
exist because agent-written documentation drifts toward two failure modes: a
feature list of sentence fragments in place of an overview, and a tone that
tells the reader what to find significant instead of stating what the code does.
The skill fixes a document shape (overview paragraph, then bullets naming real
functions and files, then a scope statement) and names the registers to avoid.
It carries no helper code and no thresholds; it is loaded for its rules and read
by whichever agent is about to write prose.

`SKILL.md` is the whole skill. Its sections are the checkable rules:

- **Structure** — open with three to five sentences of plain prose stating what
  the thing does and the problem it addresses, then bullets that name the actual
  callable or filename with one line each.
- **Tone** — descriptive and professional; no editorialising, no narration of
  the development process, no unquantified marketing adjectives, and a number
  in place of an emphatic adjective wherever one exists.
- **Omit** — bug changelogs, arrow-and-fragment shorthand in prose, and tables
  that restate an obvious directory layout.
- **Numbers and claims** — every quantitative claim must be reproducible from a
  saved artifact, and every reported threshold must say what calibrated it and
  how to re-derive it.
- **Scope statement** — close with what the thing does not do, naming the tools
  that handle it.

## Relationship to coding-standards

`coding-standards` requires a `README.md` beside every `SKILL.md` but does not
specify how to write one; it defers to this skill. The two are loaded together
when authoring a skill: `coding-standards` decides which files must exist and
where, and `doc-style` governs the prose inside them.

## Scope

Prose only. The skill does not cover code formatting, naming or line length —
those follow the conventions already present in the file being edited — and it
does not specify where files live or what a skill directory must contain, which
is `coding-standards`. Figure captions and axis labelling are governed by
`figure-style`; skill packaging, frontmatter and the `kernel.py` sidecar rules
by `skill-creator`.
