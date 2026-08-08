---
name: code-review
description: Perform thorough, actionable code reviews on files, diffs, pull requests, or code snippets. Use this skill whenever the user asks you to review code, audit code quality, check for bugs, find security issues, review a PR or diff, critique an implementation, or give feedback on code they've written or are about to merge. Also trigger when the user says things like "what do you think of this code", "anything wrong with this", "can you look over my code", "review before I merge", or "is this production-ready". Works for any language.
---

# Code Review

Perform structured, actionable code reviews that a developer can immediately act on. The goal is to be the reviewer every developer wishes they had: thorough but not pedantic, opinionated but fair, and always focused on what actually matters for the codebase.

## Philosophy

A great code review does three things: catches problems the author missed, teaches something useful, and respects the author's time. Most of what you flag should be things that would cause real pain later — bugs, security holes, maintainability traps. Stylistic nits are fine in small doses but should never dominate the review.

Think about what the code is *trying to do* before critiquing *how* it does it. If the approach itself is flawed, say so early — don't waste the author's time with line-level fixes on code that needs rethinking.

## Review process

### 1. Understand context first

Before commenting on anything, answer these questions for yourself:

- What is this code supposed to do? (Read any PR description, commit messages, docstrings, or surrounding code)
- What part of the system does it touch? (Data layer, API, UI, infrastructure, tests)
- Is this a new feature, a refactor, a bugfix, or a prototype?
- What language, framework, and conventions is the project using?

This context changes what matters. A prototype gets different feedback than production code. A performance-critical path gets scrutinized differently than a one-off migration script.

### 2. Scan for high-impact issues first

Do a quick pass looking only for things that would block a merge or cause real damage:

- **Bugs**: Logic errors, off-by-one, null/undefined access, race conditions, unhandled error paths
- **Security**: Injection vectors (SQL, XSS, command), auth/authz gaps, secret exposure, unsafe deserialization, path traversal
- **Data integrity**: Missing transactions, partial writes, unchecked inputs that hit a database
- **Breaking changes**: API contract violations, schema changes without migration, removed public interfaces

If you find any of these, flag them immediately — they're the reason the review exists.

### 3. Evaluate design and structure

Once you're satisfied there are no showstoppers, look at the shape of the code:

- Does the abstraction make sense, or is it fighting the domain?
- Are responsibilities clearly separated, or is one function/class doing too much?
- Are there obvious DRY violations (copy-paste with slight variations)?
- Would a new team member understand this code in 6 months?
- Are there missing edge cases or error handling gaps?
- Is the code testable? If it's hard to test, that's usually a design smell.

### 4. Check operational concerns

- **Performance**: Unnecessary loops, N+1 queries, unbounded memory, blocking calls in async paths, missing indexes
- **Observability**: Would you know if this code failed in production? Are errors logged with enough context?
- **Configuration**: Hardcoded values that should be configurable, environment-specific assumptions
- **Dependencies**: New dependencies that are heavy, unmaintained, or duplicating existing functionality

### 5. Review tests (if present)

- Do tests cover the happy path AND meaningful failure modes?
- Are tests actually asserting the right things (not just "it didn't crash")?
- Are test names descriptive enough that a failure message tells you what broke?
- Any flaky patterns (sleep-based timing, order-dependent tests, shared mutable state)?

If tests are missing and the change is non-trivial, note it — but suggest *which specific tests* would add the most value rather than just saying "add tests."

### 6. Note style and readability (lightly)

Only flag style issues that genuinely hurt readability or violate the project's established conventions. Do not bikeshed on formatting if the project has a formatter. Examples of worthwhile style comments:

- Misleading variable/function names
- Comments that contradict the code
- Deeply nested logic that could be flattened
- Magic numbers without explanation

## How to deliver feedback

### Severity levels

Tag every finding with a severity so the author knows what to prioritize:

- **🔴 Must fix**: Bugs, security issues, data corruption risks, broken functionality. The PR should not merge with these.
- **🟡 Should fix**: Design problems, missing error handling, performance issues, missing tests for critical paths. Strong recommendation to address before merge.
- **🟢 Consider**: Readability improvements, minor refactors, style suggestions, alternative approaches. Nice to have, author's call.
- **💬 Note**: Observations, questions, or context for future work. Not asking for changes.

### Format each finding clearly

For each issue, provide:

1. **What** the problem is (be specific — point to the exact line or pattern)
2. **Why** it matters (what could go wrong, or what makes it hard to maintain)
3. **How** to fix it (suggest a concrete alternative — code snippet when helpful)

Bad: "This function is too long."
Good: "This function handles parsing, validation, and persistence in one 80-line block. If the validation rules change (likely, given the TODO on line 42), you'd have to re-test the persistence logic too. Consider extracting `validate_input()` and `persist_record()` — each becomes independently testable."

### Tone

- Be direct but not harsh. "This will crash if `user` is None" is better than "You forgot to handle None."
- Phrase suggestions as improvements to the code, not criticisms of the author.
- Acknowledge good decisions when you see them — it reinforces what to keep doing and shows you actually read the code.
- If you're unsure about something, say so. "I'm not sure this handles the concurrent case — worth checking" is more useful than a false-positive flag.

## Output structure

Organize your review like this:

```
## Summary
[1-3 sentence overview: what the code does, overall impression, and the most important finding]

## Critical issues
[🔴 items — if none, skip this section entirely]

## Recommendations
[🟡 items]

## Suggestions
[🟢 and 💬 items]

## What's done well
[Briefly note 1-2 things the code does right — skip if nothing stands out]
```

If the code is clean and well-written, say so and keep it short. Not every review needs to be long. A review that says "This looks solid — clean separation, good error handling, tests cover the main paths. One minor thought: [suggestion]" is a perfectly good review.

## Language-specific things to watch for

You don't need to memorize these — just keep them in mind when reviewing code in these languages:

- **Python**: Mutable default arguments, bare `except:`, missing `__init__.py` in packages, type hint consistency, async/sync mixing
- **JavaScript/TypeScript**: Unhandled promise rejections, `==` vs `===`, missing `await`, stale closures in React hooks, any-typing that defeats TypeScript's purpose
- **Java**: Resource leaks (unclosed streams/connections), checked exception swallowing, mutable objects in shared state
- **Go**: Ignored errors (`_`), goroutine leaks, nil pointer on interface values
- **Rust**: Unnecessary `.clone()`, `.unwrap()` in library code, lifetime issues papered over with `'static`
- **SQL**: Injection via string concatenation, missing indexes on WHERE/JOIN columns, SELECT * in production code

## Adapting to what the user gives you

- **Full file or snippet**: Review as-is. If context is missing, note assumptions.
- **Diff or PR**: Focus on changed lines but flag if a change breaks something in surrounding code.
- **"Is this production-ready?"**: Apply all checks above with extra weight on security, error handling, and observability.
- **"Quick look"**: Prioritize only 🔴 and 🟡 items. Skip 🟢 unless something really jumps out.
- **Specific concern** ("is this thread-safe?"): Focus the review on that concern, but still flag any 🔴 issues you spot along the way.
