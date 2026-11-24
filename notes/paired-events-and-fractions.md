# Which local collections can test donor replacement and size-fraction enrichment?

Metadata inventory for P1 and P4, requested in relay message 248917669285642240. Counts below refer to field-event keys or library records as specified. No abundance outcomes were inspected.

| Collection | What is linked | What can be tested |
|---|---|---|
| CMO winters | 97 date-concordant, QC-screened DNA–RNA library pairs over 39 dates. One treated and one control pool per winter. | P1 descriptive trajectories and RNA/DNA processing checks. The number of independent, comparable donor-replacement events remains unestablished. |
| ASGARD 2017 | 31 process-event keys with 0.2, 3 and 20 µm 16S records; two further process keys with 3 and 20 µm only. | A replicated fraction contrast for candidate taxon membership. MAG-based dependence requires a separate shotgun-to-event join. |
| Amundsen 2023 | At HS22-003, one matched surface collection with three nominal pore sizes and six nonempty shotgun FASTQs. A second small/large candidate has an unresolved shotgun ID. | A concrete P4 workflow pilot. The field protocol sets the actual retained size interval. Independent events are needed for a general enrichment claim. |

## P1. Paired DNA and RNA

The dedicated raw RNA root `/matika/seqs/rnaseq/raw` contains two CMO sequencing directories. The RNA manifest's `metagenome_genohub_id` matched 112 of 130 rows to the 1,906-row metagenome manifest; 104 pairs retained blank QC or a corrected swap on both sides, and 97 of those had equal, nonblank collection dates. Details and checksums are in [the feasibility audit](vision-codex-plexi.md). Among the 104, 90 pairs are marked EXP or CTL. The other 14 consist of two river/estuary site samples on the same December 2024 date, four other samples on one December 2024 date (two MESO, two labelled SF), a kelp sample, and seven undated virus-preparation samples. These additional labels do not establish a validation series with exchange-compatible donor replacements. Tubes, fractions and dates drawn from the same pool share an exposure history.

The screen covers the dedicated RNA root and its curated CMO manifest. Older survey folders, external datasets and unpublished field sheets can contain further assays. Claiming a second independent cohort will require matching DNA/RNA to biological event and date, distinguishing independent vessels or field collections, then checking recipient coverage, donor turnover and overlapping guild capacity under a frozen protocol.

## P4. Amundsen shotgun candidates

The Micros­cape AMUN collection table has 190 rows. Its `pore.size` and `seq.filt` fields distinguish 118 nominal 0.22 µm rows coded 0 from three 0.22 µm rows coded 1, 60 rows at 3 µm and nine at 20 µm (all coded 1). The separate `prefilter.size` field matters: 350 µm mesh is upstream of many samples. Repurposing a prefilter size as the analytical fraction would produce incorrect contrasts. Match by cruise, date, time, station and target depth. Of the three 0.22 µm rows coded 1, two had a corresponding 3 µm row under this strict key.

- **HS22-003, 2023-10-16 at 16:30, surface.** The 0.22, 3 and 20 µm rows map to Genohub 9098034 plate 8 F05, plate 9 C06 and plate 9 E05, respectively. The DNA manifest reports the same collection date, station and 2.2 m depth for all three, with blank QC flags. All six R1/R2 paths exist as nonempty files under `/matika/seqs/metagenome/illumina/reads/raw/HNHG2DSXC/GC3F-UG-9920/`. Establish filter order and shared bottle ID before interpreting retained intervals. Even if collected serially, the lower fraction's upper cutoff would be determined by the preceding filter; 0.8 µm is not documented here.
- **AMD23-29, 2023-10-06 at 01:30, surface.** The 0.22 and 3 µm rows coincide. The 0.22 µm DNA ID resolves; the 3 µm sample ID did not resolve using punctuation-insensitive case-normalized equality of the AMUN `genohob.id` and DNA `sample_name`. Exclude from confirmed shotgun pairs until a primary crosswalk resolves the second library.
- **AMD23-6, 2023-10-13 at 12:30, surface.** The third candidate small fraction had no 3 µm row with the same strict key. The other 118 rows at 0.22 µm carry a different filtration code. Combining them with the three would change the contrast.

Primary tables: `/matika/microscape/out_dada_20241008/analysis_AMUN/meta.amun.complete.csv` (SHA-256 0023440f9404) and `/data/dev/genice-sample-metadata/outputs/genohub_metagenome_metadata.tsv` (ca36581041af). The first file uses the literal field spelling `genohob.id`.

## P4. ASGARD amplicon series

The ASGARD SRA export contains 313 records, all labelled AMPLICON and 16S. Of these, 247 have a 0.2 µm filter, 33 a 3 µm filter and 33 a 20 µm filter. Grouping by collection timestamp, coordinates, depth, station and Process/Survey designation gives 33 process keys and 189 survey keys. Of the 33 process keys, 31 have all three sizes and two have only the larger two. All survey keys have 0.2 µm only. Multiple libraries may share a key. Inference should block stations/casts; a depth series from one cast is not independent replication of cruise conditions. Confirm serial filtration and field replicates from the collection protocol before prescribing the size contrast.

The 16S series can test whether a prespecified membership signal transfers across fractions and independent stations. P4's dependence-per-auxotrophy outcome needs matched genome-resolved coverage. A directory labelled ASGARD-NovaSeq exists under the raw Illumina metagenome tree; the directory label or platform alone cannot establish assay identity or link those files to the individual process collections. Primary file: `/matika/microscape/out_dada_20241008/analysis_ASGARD/ASGARD2017_SRAmetadata_20260206.csv` (f3f6663b8b7f).

## Search boundaries

I screened delimited analysis tables directly under `/matika/microscape/out_dada_20241008/analysis_*` for fraction and pore-size headers and inspected both the Amundsen and ASGARD event tables. The TCN Genohub sheet has six 0.7 µm records among 96 rows and 90 blank size records. It gives no second numeric size for the proposed contrast. In the curated shotgun manifest, the Watershed project has both 0.2 and 0.7 µm records; those cutoffs alone do not document a coarse particle fraction. The raw TARA directory inspected contains accession/download lists without resident FASTQs at its top level. This is a survey of identifiable local records, not an exhaustive assay inventory for every project or a download of external data.

## External paired aquatic cohorts for P1

Salazar et al.’s companion portal describes 129 size-matched metagenome/metatranscriptome pairs in Table W6. This is a candidate for geographic replication. Check station clustering, genome assignments and focal guild coverage before counting eligible events. The 129 is the authors’ pair count. [Tara Oceans companion tables](https://ocean-microbiome.org/v2/)

Gilbert et al.’s main text reports eight DNA time points at Western English Channel station L4 across three seasons in 2008. Seven have RNA because the August 28 morning RNA sample was lost in transit. The shared site and day/night sampling make this a small temporal sensitivity series. [Gilbert et al., main text](https://doi.org/10.4056/sigs.1202536)
