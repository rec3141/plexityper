# P6 field data gate (protocol 9a47534487ba, "Data gate"): metadata only, no RNA counts opened

Checked 2026-11-20 by claude-plexi-2.

## Sources
- `/data/dev/genice-sample-metadata/outputs/cmo_rna_metadata.tsv`: the only local paired-RNA
  metadata found; 130 rows (CMO 2024-25: 66, CMO 2025-26: 64).
- `/data/dev/genice-sample-metadata/outputs/genohub_metagenome_metadata.tsv`: DNA only.
- `LWFC_2025_RosetteChemistry_Compiled_ for_ODV.csv`: 16 rows of rosette chemistry, with no
  RNA library linked.

## Result: Inconclusive. The gate is not met in the inspected local metadata.
- **Chemistry, CMO RNA.** Filled values out of 130 rows:
  - water_temp_c 40; salinity 9; pH 0; do_pct 0.
  - No NO₃⁻, NO₂⁻, NH₄⁺, SO₄²⁻, H₂S, CH₄, DIC or DOC columns.
  - The minimum chemistry in item 1 fails for every row.
- **Redox regimes.** By medium: seawater 64, sea ice 48, blank 13, sea foam 2, frost
  flowers 1, MOS culture 1, CW 1; total 130 (Astra
  249690674886524928). With no O₂ recorded, **0 blocks are documented as suboxic or anoxic**.
  This is absence of evidence, not evidence that every block was oxic. The ≥ 3 suboxic/anoxic
  block requirement is not met in the inspected metadata.
- Of the 130 RNA rows, 112 carry a metagenome_genohub_id. These IDs alone don't establish
  matched depth and time pairs, or whether a time series is feasible.
- **LWFC 2025 rosette.** It has NO₂⁻, NO₃⁻, SO₄²⁻, DIC and pH, but no O₂, NH₄⁺, H₂S or CH₄,
  and no paired RNA was found.
- **Amundsen 2023.** Only MIMARKS survey sheets were found; no paired RNA metadata turned up
  locally.

## Consequences (per protocol)
- The P6 field verdict is "Inconclusive (data gate not met)".
- External redox-stratified data with paired chemistry and RNA would need grid's choice.
- The precondition test, a joint DNA/RNA state-space model on CMO, still runs whatever the gate
  says. It needs per-MAG DNA and RNA coverage time series for CMO.
