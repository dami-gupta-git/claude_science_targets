# Custom skills

The thirteen skills written for this account, with a plain-language description
of each. Source for all of them lives in this directory, one folder per skill.

Most answer some version of the same question — is this gene or compound worth
pursuing? — and the recurring theme is that the obvious analysis gives the wrong
answer unless a specific control is run first. Three are writing, coding, and
review conventions rather than analysis tools, and one is a personal daily
routine.

Entries are grouped by what they are for. They are not ordered by date: the
files carry the timestamp of the last catalog sync rather than when they were
written, so authorship dates cannot be read off the filesystem.

## Where things live

`claude_science/` holds two directories. `skills/` is this one, with a folder
per skill containing `SKILL.md`, `README.md`, `kernel.py` and `tests/`.
`results/` holds everything an analysis produced, laid out
`results/<topic>/<run>/` with the code that produced a run in that run's
`scripts/` subfolder:

```
results/target_triage/usp1/          tables, figures and README for one run
results/target_triage/usp1/scripts/  the code that produced them
results/protein_structure/wrn/
```

Skill directories are lowercase-with-hyphens, because the Agent Skills standard
requires the directory name to match the `name:` field in `SKILL.md`
frontmatter. Result directories are `snake_case`, so the separator alone says
which tree a path belongs to. Results are never written into a skill directory:
publishing uploads the whole folder, so a stray table or figure would ship to
the skill registry and appear on every machine that syncs the catalog.

These skills are also published to the Claude Science registry, which is what
makes them loadable in a session. This directory is the source you edit; the
registry copy is generated from it.

## Judging a drug target

Four of these form a stack. `target-triage-public-data` sets the order of work,
`depmap-local` and `opentargets-evidence` supply one half of the evidence each,
and `depmap-fusion` joins them into a verdict.

### target-triage-public-data

The end-to-end workflow for judging whether a gene is a plausible drug target,
deliberately ordered so the cheapest evidence that could rule it out runs first —
a gene that cancer cells tolerate losing is not a single-agent target, and no
amount of later analysis changes that. It runs from genetic dependency through
druggability, then drug-response correlations, then patient-survival data. Its
main contribution is naming the trap in the drug-response step: genes tied to how
fast cells divide appear to correlate with sensitivity to *every* cytotoxic drug,
so an apparently drug-specific finding is often nothing of the kind. The skill
corrects for this and reports the corrected number.

### depmap-local

Answers "do cancer cells actually need this gene?" from a large public experiment
held on disk. The underlying data covers roughly 1,200 cancer cell lines in which
each of about 18,000 genes was switched off in turn, and this skill reads it
without loading the whole thing into memory. It reports whether a gene is needed
by nearly every cell line — real, but a poor drug target, since healthy tissue
needs it too — or only by a specific subset, and which cancer types or genetic
backgrounds those are. It also handles the file-format traps that silently return
empty results, and covers the companion screen of drug responses across the same
cell lines.

### opentargets-evidence

Pulls human evidence for a gene–disease link from the Open Targets database: how
strongly the two are associated, what kind of evidence supports it (genetic
studies carry far more weight than counting how often the two are mentioned
together in papers), whether the protein is the sort of thing a drug can
realistically bind, and whether any drugs against it have reached the clinic. It
documents the specific query shapes that work against the live API, which is the
part that costs time to rediscover.

### depmap-fusion

Combines the two evidence sources above and reports where they disagree. Human
evidence says nothing about direction — so tumour-suppressor genes, where
blocking them would make things worse, rank near the top of cancer lists. The
cell-line data has direction but no human evidence. Fusing them sorts each
candidate into one of ten verdicts, from "both sources agree, strongest case"
through "switching this gene off actually helps the cancer grow". Across two
cancers tested, about two-thirds of the top-ranked candidates turned out not to
be viable targets.

### marker-contrast-null

Checks whether a "this drug target only matters in patients with mutation X"
claim is real or a statistical accident. The trap it guards against: if you pick
one mutation out of thousands and test it, you will often find an
impressive-looking result purely by chance. This skill re-runs the same test
against every other candidate mutation and reports where your chosen one ranks —
so a result in the top 4% of 1,700 candidates is correctly read as unremarkable
rather than as a discovery. It also checks two ways the result could be an
artefact: whether the mutant cell lines are simply frailer overall, and whether
unrelated genes show the same pattern.

### depmap-genetea

Answers two questions the other cell-line skills do not: which other genes behave
like this one across the panel, and what a list of genes has in common when
stated in words. The first is guilt-by-association — if switching off gene A and
gene B harms the same cell lines and spares the same ones, they likely work
together — which proposes protein complex members and pathway partners. The
second takes an unnamed set of genes, such as the hits from a screen, and names
the theme they share, tracing each term back to the reference sentence that
produced it rather than returning an unexplained category label.

## Molecules and structures

### get-protein-structure

Fetches a usable 3D protein structure starting from whatever identifier you
have — a gene name, a database accession, a protein's common name. Rather than
returning the first match, it ranks all available experimental structures by
resolution, how much of the protein they actually cover, and whether a drug-like
molecule is bound, then writes a cleaned-up file ready for downstream docking or
design work. It also reports which parts of the protein are missing from the
structure, flagging gaps that sit inside the pocket you care about. When no
experimental structure exists it falls back to a predicted one with its
confidence scores attached.

### boltz-affinity-triage

Ranks candidate drug molecules against a protein target using AI structure
prediction, with the controls needed to make the ranking mean something. The
problem it solves: these models return a confident-looking binding-strength
number for any molecule you give them, and on its own that number cannot be
trusted — partly because the model was trained on most of the public data you
might check it against. The skill scores known active compounds and
carefully-matched inactive ones alongside your candidates, so you can see whether
the model is separating real binders from decoys on this target at all.

## Evidence and claims

### ai-origination-audit

Separates what a company claims about AI's role in discovering a drug from what
the public record actually shows. Rather than recording "AI-designed" as a
property of a drug candidate, it records which specific step AI was used for —
finding the target, generating the molecule, optimising it, designing the trial —
along with the source of that claim and how strong that source is. The output is
a fixed set of tables, so candidates audited months apart stay comparable and can
be merged.

## Conventions

### coding-standards

House rules for writing analysis code. The central point is that code which
exists only inside a running session is a result nobody can reproduce — anything
worth running twice gets written to a file and saved. It sets the repository
layout described above: skill source authored into this directory first and then
published, results under `results/<topic>/<run>/` with code in `scripts/`, and
never results inside a skill folder. It covers what makes a test meaningful (a
test that still passes when you deliberately break the code is documentation,
not a check), and states that git commands are never run unless you ask.

### doc-style

House rules for writing prose — READMEs, reports, documentation. It asks for a
plain-language opening paragraph rather than a feature list, a descriptive tone
that reports what the data show and leaves the reader to judge importance, and no
marketing adjectives. It also asks documents to describe how something currently
behaves rather than narrating its repair history, and requires every number in a
document to be traceable to the file that produced it.

### code-review

Structured, actionable code review of a file, a diff, or a pull request. The
priority is problems that cause real pain later — bugs, security holes,
maintenance traps — over style preferences, and each point is meant to be
specific enough to act on immediately. Works for any language.

## Personal

### morning

Renders a morning brief as a single styled page: the shape of the day, what needs
attention, and what is already handled. It can also be set up as a recurring
weekday task. It runs only on explicit request — a passing question about the
day's schedule is answered directly rather than triggering the full brief.
