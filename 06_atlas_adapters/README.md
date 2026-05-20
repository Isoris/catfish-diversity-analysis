# 06_atlas_adapters

Adapter scripts that transform pipeline outputs into the four
schema-conforming JSON payloads consumed by the **diversity-atlas**:

| Atlas page | JSON file (in diversity-atlas/data/)   | Adapter                            | Upstream pipeline                                         |
| ---------- | -------------------------------------- | ---------------------------------- | --------------------------------------------------------- |
| `texture`  | `texture_metrics.json`                 | `adapt_texture_metrics.py`         | [04_window_diversity_texture/](../04_window_diversity_texture/) (ready) |
| `roh`      | `roh_gene_overlap.json` (extension)    | `adapt_roh_gene_overlap.py`        | [05_roh_metrics_and_genes/](../05_roh_metrics_and_genes/) (ready)       |
| `divergence` | `divergence_network.json`            | `adapt_divergence_network.py`      | pairwise FST / DXY pipeline (**not yet shipped**)         |
| `burden`   | `functional_burden.json`               | `adapt_functional_burden.py`       | bcftools csq + snpEff + VESM (**not yet shipped**)        |

Schemas live in
[diversity-atlas/atlases/diversity/registries/schemas/](../../diversity-atlas/atlases/diversity/registries/schemas/):
`texture_metrics_v1.schema.json`,
`roh_gene_overlap_v1.schema.json`,
`divergence_network_v1.schema.json`,
`functional_burden_v1.schema.json`.

## Driver

Run everything in one shot:

```bash
bash STEP_C01_emit_atlas_payloads.sh
```

The driver sources `../00_config.sh` for paths, runs each adapter, and
copies the emitted JSONs into `${OUTBASE}/13_atlas_payloads/`. To
publish to the atlas, copy those files into
`diversity-atlas/data/` (or set `--copy-to-atlas` on the driver to do
it in one step).

## Two-tier adapter behaviour

Each script handles missing upstream inputs the same way:

1. **All inputs present** → emit a fully populated JSON with
   `version: "v1-derived"`, `generated_at: <ISO timestamp>`, and the
   real data.
2. **Inputs missing** → emit a schema-valid stub with
   `version: "v0-stub"`, `generated_at: <ISO timestamp>`, empty arrays,
   and a `_doc` field naming what's still pending. The atlas page
   detects the empty payload and renders a "data pending" card; nothing
   crashes.

This means the driver always emits valid JSON. Plug in the missing
upstream pipelines later without touching the atlas — the next driver
run picks up the new inputs.

## Conventions

- Adapters write to stdout if no `--out` is given. The driver always
  passes `--out`.
- Inputs are TSVs / BEDs produced by sibling modules under
  `${OUTBASE}/`. Each adapter documents its input contract in its
  `--help`.
- All scripts are Python 3 + pandas + numpy only (matches the rest of
  the repo's dependency surface).

## Adding a new payload

Mirror the existing adapters:

1. Drop the schema JSON in `diversity-atlas/atlases/diversity/registries/schemas/`.
2. Write `adapt_<name>.py` in this folder. Accept `--out` + per-input
   args. Support the two-tier behaviour above.
3. Add a line to `STEP_C01_emit_atlas_payloads.sh`.
4. Update the table above.
