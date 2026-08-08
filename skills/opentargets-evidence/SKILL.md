---
name: opentargets-evidence
description: Pull target-disease evidence from the Open Targets Platform GraphQL API - association scores by data type, tractability by modality, known drugs and max clinical stage, genetic constraint, safety liabilities, and the prioritisation panel. Also ranks targets for a disease (EFO/MONDO id). Use when asked what human evidence supports a target, whether it is druggable or has approved drugs, or which targets are top-ranked for an indication.
---

# Open Targets evidence

Reached through the **clinical-genomics connector**, from the `repl` tool only:

```python
host.mcp("clinical-genomics", "open_targets_graphql", query=Q, variables={...})
```

Convenience wrappers exist for two common cases - `open_targets_disease_targets`
(ranked targets for a disease) and `open_targets_drug` (drug by ChEMBL id).
Everything else goes through `open_targets_graphql`.

## Schema shapes verified against the live API

These bite if guessed - all confirmed working:

- `tractability` fields are **`label` / `modality` / `value`**. There is **no `id`**.
  Modality codes: `SM` small molecule, `AB` antibody, `PR` PROTAC/degrader,
  `OC` other clinical. Filter on `value == true`; the raw list includes false rows.
- `drugAndClinicalCandidates` takes **no arguments** - do not pass `size`.
  It replaced the removed `knownDrugs` field.
- `associatedDiseases` paginates with `page:{index,size}`.
- `geneticConstraint` returns rows per `constraintType` (`syn`/`mis`/`lof`);
  the `lof` row carries `upperRank` (lower rank = more constrained).
- `safetyLiabilities` is often an empty list - absence is not safety.

## Resolving a symbol to an Ensembl id

`target(ensemblId:)` needs an ENSG. Resolve first:

```graphql
query($q:String!){ search(queryString:$q, entityNames:["target"], page:{index:0,size:5}){
  hits{ id object{ ... on Target { approvedSymbol biotype } } } } }
```

Take the hit whose `approvedSymbol` matches exactly and whose biotype is
`protein_coding` - a plain search for `KRAS` also returns the pseudogene `KRASP1`.

## Dossier query

One round trip for the full target picture:

```graphql
query($id:String!, $n:Int!){
  target(ensemblId:$id){
    id approvedSymbol approvedName biotype
    tractability{ label modality value }
    safetyLiabilities{ event datasource }
    geneticConstraint{ constraintType score upperRank oe }
    associatedDiseases(page:{index:0,size:$n}){
      count rows{ score disease{ id name therapeuticAreas{ id name } }
                  datatypeScores{ id score } } }
    drugAndClinicalCandidates{ count rows{ maxClinicalStage drug{ id name drugType } } }
    prioritisation{ items{ key value } }
  }
}
```

## Interpreting scores

- The overall association score is an evidence-weighted aggregate in 0-1;
  read `datatypeScores` to see **why** (`genetic_association` and
  `known_drug` are far stronger prior evidence than `literature`, which is
  co-mention-driven and inflates for well-studied genes).
- A high score means "human evidence links this gene to this disease" - **not**
  that inhibiting it helps. Direction is not encoded: tumour suppressors such
  as TP53 rank at the top of cancer lists. Pair with a functional KB before
  calling anything a drug target (see `depmap-fusion`).
- `prioritisation` returns signed values in [-1, 1]; the number of keys is PER-TARGET, not fixed (KRAS returns 17, WRN 15) - iterate the items, never index by position
  (constraint, pockets, ligands, membrane/secreted, mouse KO, paralogs).

## Rate limits

The connector is shared across subagents. Batch fields into one dossier query
rather than looping; a 30-target sweep at one call each takes ~25 s.
