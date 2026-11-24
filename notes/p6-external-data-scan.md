# P6 external data scan (2026-11-24, claude-plexi-2; grid asked "look for some data that will work for P6. tara oceans maybe? or GEOTRACES")

Gate (P6 9a47534487ba plus table notes §1, §5.3):
- paired DNA + RNA from the same water;
- minimum chemistry: O2, NO3, NO2, TAN, SO4, ΣH2S, CH4, DIC, DOC, T, pH;
- ≥ 20 events, ≥ 6 independent casts, ≥ 3 oxic and ≥ 3 suboxic/anoxic blocks.

Sources were read as web summaries and abstracts only. No data were downloaded and no RNA was opened.

| Candidate | Paired metaG + metaT | Independent blocks | Chemistry vs minimum | Verdict |
|---|---|---|---|---|
| **Saanich Inlet** (Hawley 2017 Sci Data 4:170160; Torres-Beltrán 2017 4:170159) | 90 metaG, 62 metaT (46 unique), 6 depths 10–200 m, across monthly cruises 2006–2014; pairing per sample in Supp. Table 2 (not yet read) | Likely ≥ 6 cruises; each spans the oxycline (oxic 10 m, suboxic/anoxic/sulfidic deep) | Has O2, NO3, NO2, NH4, H2S, CH4, N2O (tier open until units and paired-event coverage are checked), CO2, T, S. **Missing DIC, DOC and pH**; SO4 can be derived from salinity | **Eligibility unverified; bottle/cast and chemistry audit pending.** Strongest candidate |
| Cariaco Basin MetaOmics (BCO-DMO awards 562424/562429) | 48 metaG + 48 metaT (free-living / particle fractions), 6 depths | **2 cruises only**; size fractions are not independent blocks | CARIACO time series has O2, NO3, NO2, NH4, H2S (DIC, pH and CH4 unconfirmed) | Fails the ≥ 6-block gate. Could join a pooled cohort |
| Tara Oceans (Salazar 2019 Cell; PANGAEA env. context) | 187 metaT / 370 metaG globally; mostly surface/DCM; few mesopelagic OMZ metaT | One cast per station; few OMZ stations with RNA | NO2, NO3+NO2, PO4, Si, O2, carbonate (CT/AT, matched within 100 km / 2 days); **no NH4, H2S, CH4; DOC Arctic only** | Fails chemistry and anoxic blocks |
| GEOTRACES / bioGEOTRACES (Biller 2018 Sci Data) | **480 metaG, no metaT found** | many | Excellent chemistry | Fails: no RNA. Useful for DNA-only work (P2/P3) |
| ETNP/ETSP OMZ (Stewart 2012; Glass 2015; Ganesh 2015) | ≤ 9 metaT, 1–3 stations | < 6 | NO3, NO2, O2; little NH4, H2S, CH4 | Fails |
| Lake Erie (Brown 2025 MRA) | 2 + 2 | 1 | partial; hypoxic, not anoxic | Fails |
| Lake Yojoa, Honduras (Zenodo 8050258) | metaT genome-resolved; 3 sites × 2 dates × 2 depths | ≤ 6; anoxic only in June | unverified | Unverified second candidate |

## Recommendation
Audit Saanich Inlet first; eligibility is unverified. If the audit passes, it would need P6-A1, a chemistry addendum fixed before any RNA is opened, plus cross-family review and grid's go:
- **DIC/HCO3⁻:** from measured dissolved CO2 plus a pH source (DFO or regional data, or a stated sensitivity range).
- **DOC:** a fixed regional value with sensitivity range. Correction (Astra 251077299830964224): after normalisation, a fixed DOC also shifts the modelled heterotroph shares relative to the other pathways. So it is a substantive assumption, to be stated as such in the addendum and given a sensitivity range.
- **SO4²⁻:** from salinity, using the conservative seawater ratio.

Next steps:
1. Read Supp. Table 2 (the cruise × depth pairing). This is a small xlsx download, so it needs grid's OK.
2. Check the Dryad chemistry (doi:10.5061/dryad.nh035) for the same cruises.
3. Data access: JGI IMG study "Marine microbial communities from expanding oxygen minimum zones in the northeastern subarctic Pacific Ocean"; SRA SRP043213.

## Audit items before any eligibility claim (Astra 251077299830964224)
Saanich eligibility under frozen P6 is **unverified** until a metadata audit shows each of the following:
- same-bottle DNA/RNA pairing. The 46 unique RNA samples plus 16 replicates have to be collapsed to events.
- compatible size fractions.
- distinct cast IDs. Torres-Beltrán used the deepest cast's downcast CTD profile, so cruise × nominal depth cannot certify cast provenance.
- complete chemistry. The Dryad dictionary lists headspace CO2 and no pH, DIC or DOC.
- CH4 and N2O units, checked in the actual file before assigning chemistry tiers, because the published labels conflict (nM vs µM).

## Further audit and addendum items (Carbon 251077656854315008)
1. **SUP05 nitrate-coupled sulfide oxidation.** In the suboxic and sulfidic zone, SUP05 couples HS⁻ oxidation to NO₃⁻ reduction.
   - Their nap/nar/nir/nor transcripts would land on DEN_* rows, which model acetate as the donor.
   - Pre-register the frozen SOX_NO3 MAG-resolved split alongside SOX_AER, and state the donor-mismatch caveat (table notes §2.1).
2. **pH and carbonate must vary with depth.** Carbon's figure of pH ~8.0 → ~7.4 from 10 to 200 m is an **illustrative sensitivity scenario, not a measured S3 profile**. The audited archive has no pH, and no dated, depth-resolved source or pH scale is cited (Astra 251077834248208384).
   - A 0.6-unit difference would shift HCO₃⁻/CO₂ about 4×. That explains why pH matters, but it does not establish the profile for the paired casts.
   - Look for DFO IOS Saanich / Line-P carbonate casts first.
   - Otherwise reconstruct pH(z) and HCO₃⁻(z) from Henry's-law CO₂(aq) plus TA(S), corrected for HS⁻ and NH₃, with a pH sensitivity.
   - DOC also changes ratios *within* the heterotroph rows: SRB_ORG has K_s(acetate) = 5.0 µM vs 15.6 µM for AER_HET/DEN/DNRA.
3. **Temporal dependence at one station (S3).** Deep water is renewed roughly once a year, so cruises within one stagnant period share water masses.
   - Add a leave-one-renewal-year-out check alongside leave-one-cruise-out.
   - Check that metaT library prep and rRNA depletion were consistent across 2006–2014.
4. **Backup candidate: Baltic Landsort or Gotland Deep redoxcline.** Carbon's claim that the IOW monitoring chemistry includes O₂, nutrients, H₂S, CH₄, pH, DIC and DOC, together with multi-cast metaG+metaT, is **unverified** (Opus has not checked it).

## Saanich Supp. Table S2 read (2026-11-24; grid saved the files; inventory only, no RNA opened)
Files are in data/external/saanich/download_20261124/hawley_supp/: MOESM151 (Table S2) 732686f9cf36, MOESM150 (Table S1) f668abb12ebc, and the ISA zip MOESM149 4b8be87b9efb.
- Table S2 has 119 rows: 90 metaG, 46 metaT. **44 sample IDs carry both** (the same SampleID, i.e. cruise × depth). Same-bottle pairing is not yet established.
- All 44 are from station SI03.
- **8 cruises** (2010–2012) have both. Each cruise is one potential block, which is ≥ 6:

  | Cruise | Date | Depths (m) |
  |---|---|---|
  | 42 | 2010-Feb | 10, 100, 120, 135, 150, 200 |
  | 47 | 2010-Jul | 100, 120, 135, 150, 200 |
  | 48 | 2010-Aug | 10, 100, 120, 135, 150, 200 |
  | 54 | 2011-Feb | 100, 120, 135, 150, 200 |
  | 72 | 2012-Aug-1 | 10, 100, 135, 150, 165, 200 |
  | 73 | 2012-Aug-28 | 10, 165, 200 |
  | 74 | 2012-Sep-10 | 10, 100, 120, 135, 150, 165, 200 |
  | 75 | 2012-Sep-20 | 10, 100, 120, 135 (×2), 200 |

- **3 calendar years (2010–2012).** These are not verified renewal-cycle blocks; the renewal timing still has to come from the O2/density record (Astra bb339d5d1ba1).
- **Still open:**
  - Oxic and anoxic classification per block needs O2 from the Dryad chemistry.
  - Cast IDs.
  - Size fractions: check the ISA assay files.
  - pH, DIC and DOC are absent.
  - CH4/N2O units.

## Saanich chemistry join (Dryad zip saved by grid; 899e821a2bb4, files unzipped and hashed in place)
- The README gives **CH4 and N2O in nanomolar**, which resolves the unit conflict.
- There are 43 unique paired cruise × depth events (44 sample IDs; cruise 75 at 135 m appears twice).
- **Block gate: PASS on metadata.** 8 cruises, each with at least one oxic event (CTD_O2 ≥ 10 µM) and at least one suboxic/anoxic event. That gives 8 oxic and 8 anoxic blocks against a requirement of 3 each, and 43 events against 20.
- **Chemistry, out of 43 events:**
  - Complete: CTD_O2, NO3, NH4 and NO2 (43/43).
  - Nearly complete: H2S, CH4 and N2O (42/43; the missing event is cruise 73 at 200 m).
  - **No headspace CO2 for any paired event.** pH, DIC and DOC are absent throughout.
  - SO4 can come from CTD salinity.
- **Not yet checked:**
  - same-bottle DNA/RNA pairing;
  - cast IDs (the chemistry uses the deepest cast);
  - size fractions: per Hawley's Methods, metaT used 0.22 µm Sterivex with an in-line 2.7 µm GF/D prefilter; metaG has still to be confirmed as the same;
  - the repeated cruise 75 / 135 m entry (two sample IDs, one event?).
- **Consequence:** P6 on Saanich needs P6-A1 to state how DIC/HCO3⁻, pH and DOC are supplied. Measured CO2 is not available, so this would be an assumption-driven sensitivity analysis rather than measured chemistry.

## Provenance audit result (P6-A1 §1, 2026-11-24)
See data/external/saanich/biosample_20261124/AUDIT.md (local only).
- Table S2's metaG column is misaligned for cruises 72–74, so pairing was rebuilt from
  BioSample titles.
- Proxy roster: 43 events, 8 blocks (8 oxic, 8 anoxic).
- Strict roster: 0 events; bottle and cast are undocumented.
- Fractions are compatible (0.22–2.7 µm).
- The cruise 42 DNA assignment is inferred from 'F_10'.

## MGnify check (grid's suggestion, 2026-11-24)
- **MGYS00005865** ("EMG produced TPA metagenomics assembly of PRJNA247822", pipeline v5.0) has
  141 assembly analyses over 138 samples.
  - It covers **all 38 roster DNA BioSamples**, and 49 of 54 RNA.
  - Each analysis offers: processed contigs, predicted CDS (aa/nt), a combined GFF
    (eggNOG/InterPro/antiSMASH), a per-CDS eggNOG TSV with KEGG KOs, a KEGG orthologues
    table, InterPro and Pfam.
- MGYS00005850 is the raw SRP043213 study; its v5 analyses are amplicon.
- **No MAGs:** the MGnify genome catalogues have no PRJNA247822 genomes. Marine v2.0 (13,223
  species) holds cross-study representatives, not per-sample bins.
- **Consequence:** the per-event assembly step (the most expensive part of the run spec) can
  be replaced by the MGnify per-sample assemblies with their per-CDS KO annotation. Binning
  would still be needed for the binned-KO classifier, or the primary could treat all contigs as
  unbinned (unambiguous KOs only). Either way is a run-spec change that needs cross-family review.
