# code-review

Conventions for reviewing code — a file, a snippet, a diff or a pull request —
in any language. The skill orders the review so that high-impact findings are
looked for before line-level ones: context first, then bugs and security and
data-integrity risks, then design, then operational concerns, then tests, then
readability. Each finding carries a severity tag and states what the problem is,
why it matters, and a concrete alternative, so the author can act on it without
a second round of questions. It carries no helper code and no thresholds; the
skill is prose the reviewing agent follows.

`SKILL.md` is the whole skill, organised as a review process followed by
delivery conventions:

- **Review process** — six ordered passes: establish what the code is for and
  what it touches; scan for merge-blocking issues (logic errors, injection
  vectors, auth gaps, secret exposure, unchecked writes, breaking API or schema
  changes); evaluate abstraction and separation of responsibilities; check
  performance, observability, configuration and new dependencies; assess whether
  tests assert the right things and avoid flaky patterns; note only the style
  issues that hurt readability or break the project's own conventions.
- **Severity levels** — four tags, from must-fix (the PR should not merge) and
  should-fix (address before merge) through consider (author's call) to note
  (observation, no change requested).
- **Finding format** — what, why and how for each item, with the exact line or
  pattern named rather than the category.
- **Output structure** — summary, critical issues, recommendations, suggestions,
  and what the code does well; empty sections are dropped rather than filled.
- **Language-specific checks** — the recurring traps in Python, JavaScript and
  TypeScript, Java, Go, Rust and SQL, as a reminder list rather than a
  checklist to run in full.
- **Adapting to the request** — how a full-file review, a diff, a
  production-readiness question, a quick look, and a single named concern each
  change what gets reported.

Review length is set by what the code warrants. A change with no findings beyond
one suggestion is reported as such rather than padded to fill the output
sections.

## Scope

Reading and judging code, not changing it: the skill produces findings, and
applying them is left to the author. It does not run linters, formatters, type
checkers or test suites, and it does not execute the code under review, so
findings that depend on runtime behaviour are stated as things to check rather
than as confirmed defects. Where code lives, how it is structured for reuse, and
what must be true before an analysis is called finished belong to
`coding-standards`; prose conventions for the review write-up itself, if it is
saved as a document, belong to `doc-style`.
