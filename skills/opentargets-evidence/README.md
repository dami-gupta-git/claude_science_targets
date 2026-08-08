# opentargets-evidence

[Open Targets](https://platform.opentargets.org) aggregates human genetics,
clinical trials, pathway data and literature into a target-disease association
score from 0 to 1. This skill is a query layer over the live GraphQL API,
reached through the clinical-genomics connector, and it answers two questions
about a gene: what human evidence links it to a disease, and whether it is
druggable or already has drugs against it. It also ranks targets for a disease
given an EFO or MONDO id. There are no local files — every call is fetched
live.

## Functions

All calls run in the `repl` tool only.

- `host.mcp("clinical-genomics", "open_targets_graphql", query=Q, variables={...})`
  — arbitrary GraphQL against the Open Targets schema; the entry point for
  target dossiers, symbol resolution and anything without a wrapper.
- `host.mcp("clinical-genomics", "open_targets_disease_targets", efo_id=..., size=...)`
  — ranked targets for a disease id.
- `host.mcp("clinical-genomics", "open_targets_drug", ...)` — a drug record by
  ChEMBL id.

## Queries

Ranked targets for a disease:

```python
host.mcp("clinical-genomics", "open_targets_disease_targets",
         efo_id="MONDO_0005061", size=5)
# lung adenocarcinoma, 8760 associated targets:
# EGFR 0.822 | TP53 0.740 | KRAS 0.738 | BRAF 0.702 | STK11 0.696
```

`target(ensemblId:)` needs an ENSG, so resolve the symbol first and take the
`protein_coding` hit whose `approvedSymbol` matches exactly:

```python
host.mcp("clinical-genomics", "open_targets_graphql", query="""
query($q:String!){ search(queryString:$q, entityNames:["target"], page:{index:0,size:3}){
  hits{ id object{ ... on Target { approvedSymbol biotype } } } } }""",
  variables={"q":"KRAS"})
# ENSG00000133703 KRAS protein_coding
# ENSG00000220635 KRASP1 processed_pseudogene
```

Tractability and existing drugs:

```python
host.mcp("clinical-genomics", "open_targets_graphql", query="""
query($id:String!){ target(ensemblId:$id){ approvedSymbol
  tractability{label modality value} drugAndClinicalCandidates{count} } }""",
  variables={"id":"ENSG00000133703"})
# KRAS modalities: AB, PR, SM | drugs: 3
```

## Schema constraints

Verified against the live API:

- `tractability` fields are `label` / `modality` / `value`. There is no `id`.
  Modality codes: `SM` small molecule, `AB` antibody, `PR` PROTAC/degrader,
  `OC` other clinical. Filter on `value == true` — the raw list includes false
  rows.
- `drugAndClinicalCandidates` takes no arguments; do not pass `size`. It
  replaced the removed `knownDrugs` field.
- `associatedDiseases` paginates with `page:{index,size}`.
- `geneticConstraint` returns rows per type (`syn`/`mis`/`lof`); the `lof` row
  carries `upperRank` (lower = more constrained).
- `safetyLiabilities` is often an empty list — absence is not safety.

## Interpreting the score

A high score means a gene is *involved* in a disease. It encodes no direction,
so it does not mean "inhibit this". Some genes cause disease when they are too
active (oncogenes — blocking them helps); others cause disease when they stop
working (tumour suppressors — blocking them makes things worse). Both score
highly: `TP53`, a tumour suppressor, ranks #2 for lung adenocarcinoma (0.740),
and inhibiting it would be actively harmful. Direction requires functional
data, which is what `depmap-fusion` supplies.

Read `datatypeScores` to see why a target scored as it did.
`genetic_association` and `known_drug` are much stronger prior evidence than
`literature`, which is co-mention-driven and inflates for well-studied genes.

## Rate limits

The connector is shared. Batch fields into one dossier query rather than
looping — a 30-target sweep at one call each takes ~25 s.

## Scope

This skill reports human association evidence and druggability only. It does
not score directionality or cell-line dependency (`depmap-fusion`), resolve
protein structures or binding pockets (`get-protein-structure`), or retrieve
trial records beyond the counts and maximum clinical stage that Open Targets
itself carries. The full dossier query and the cross-knowledge-base workflow
live in the `depmap-fusion` README.
