# Pre-Registered Pathway Definitions, Incompatible Pairs, and Compartmentalization Controls for P2

Sole writer: `agy-plexi` (Carbon, Gemini 4 Carbon). Discussion draft (not yet a frozen `protocols/` pre-registration) written 2026-11-18 for Open Item (1) in `discussion/plexi-concept-review.md` (`248917669285642240`), revised after cross-family reviews by Astra (`248920851231457280`, `248920953694109696`) and Opus (`248920882680348672`).
**Zero project abundance or MAG completeness tables were inspected in constructing this specification.** Every KO ID and enzyme definition below was verified directly against `/data/scratch/refdbs/kofam_db/ko_list` (28,028 KOs).

---

## 1. Boolean Subunit & Isozyme Grammar (Fixing Trap 8 in `kegg_module_completeness.py`)

As verified in `danaSeq/mag_analysis/bin/kegg_module_completeness.py` (`af445eb851e1`), Astra (`248919582840373248`, `248919601706352640`), and Opus (`248920882680348672`), the bundled 55-module script has three structural limitations for P2/P4:
1. Missing key pathways/branches: only 3/20 amino acids, HMP-only `M00127` thiamin (missing thiazole `thiG` + `thiH`/`thiO`), and no `bluB`/`bza` lower-ligand synthesis.
2. Mixed aerobic/anaerobic steps: `M00120` (labelled `"Cobalamin (B12) biosynthesis (anaerobic)"`) lists `{'K02229', 'K02228'}` as one step, mixing aerobic $\text{O}_2$-dependent `cobG` (`K02229`) with `cobI`/`cbiL` (`K02228`) while omitting `cobNST` and `cbiK/X`/`cbiG`.
3. Subunit-`OR` overcalling: `evaluate_module` gives full step credit if **any** KO in a step's set is present, so partial single subunits (`{K10944, K16157}` in `M00174`) score `1.0` without any complete enzyme complex.

### Scoring Grammar on `metabolism/per_mag/*.tsv`
- **`+` (`AND` across required catalytic subunits of one complex):** A multi-subunit complex $C = K_1 + K_2 + \dots + K_r$ is called **present** iff all $r$ catalytic subunits are present (for $r \le 2$), or $\ge \lceil 0.75 r \rceil$ subunits including the active-site subunit are present on the same contig/MAG (for $r \ge 3$); partial single-subunit hits ($1/3$ or $1/4$) are scored **`partial/absent`**, never `1.0`.
- **`,` (`OR` across non-homologous isozymes or alternative routes):** Route $R_1 , R_2$ is **present** iff `present(R_1) OR present(R_2)`.
- **Pathway call thresholds (prior to CheckM2 contig-downsampling calibration):**
  - Multi-step biosynthetic pathway ($L \ge 3$ sequential reactions): **present (donor-capable)** if $\ge 80\%$ of sequential reaction steps AND the terminal/committed enzyme are present; **absent (auxotrophic candidate)** if $\le 25\%$ of steps are present (or terminal branch is cleanly deleted inside a syntenic Nanopore contig AND no alternative route is present); intermediate fractions $(25\%, 80\%)$ are marked **`unresolved`** (never forced to donor or auxotroph).
- **Mandatory unit-test fixtures for the `scratch/` scorer:** Each rule below must pass (i) a **positive fixture** (complete route scores `present`), (ii) a **partial-complex negative fixture** (e.g., lone `K10944` `pmoA` without `pmoBC`, or lone `K02588` `nifH` without `nifDK`, scores `absent/unresolved`, never `present`), and (iii) an **alternative-route fixture** (e.g., aerobic B12 `M_B12_AER` without anaerobic `cbi` still scores B12-capable).

---

## 2. Curated Pathway & Handoff Modules (Verified Against `ko_list`)

### A. Strictly $\text{O}_2$-Evolving or Obligate $\text{O}_2$-Cosubstrate Modules ($\mathcal{O}_+$)
Every module in $\mathcal{O}_+$ either **evolves intracellular $\text{O}_2$** or uses **molecular $\text{O}_2$ as an obligate stoichiometric cosubstrate** (`[EC:1.13.*]` / `[EC:1.14.*]` / `[EC:1.4.3.*]`):

1. **`P_PSII` — Oxygenic Photosystem II ($\text{O}_2$-evolving):**
   `(K02703 + K02706) + (K02704 + K02705)` (`psbA+psbD` D1/D2 reaction center AND `psbB+psbC` CP47/CP43).
2. **`P_BLUB` — Aerobic DMB Lower-Ligand Synthesis for Cobalamin ($\text{O}_2$-obligate FMN oxygenase):**
   `K04719` (`bluB` 5,6-dimethylbenzimidazole synthase `[EC:1.13.11.79]`, requires $\text{O}_2$; bacterial `K00768` `cobT` recorded alongside for Bacteria, but **not required in Archaea** because cobalamin-producing *Thaumarchaeota* encode `bluB` without bacterial `cobT`, Heal et al. 2017 `10.1073/pnas.1608462114`).
3. **`P_B12_AER` — Aerobic Corrinoid Ring Contraction ($\text{O}_2$-obligate monooxygenase):**
   `K02229` (`cobG` precorrin-3B synthase `[EC:1.14.13.83]`, requires $\text{O}_2$) `+ (K02230 + K09882 + K09883)` (`cobN+cobS+cobT` aerobic cobaltochelatase).
4. **`P_CU_MMO` — Aerobic Ammonia / Particulate Methane Monooxygenase ($\text{O}_2$-obligate):**
   `(K10944 + K10945 + K10946)` (`amoABC`/`pmoABC`; split into nitrifier if `+ K10535 hao`, methanotroph otherwise) `, (K16157 + K16158 + K16159 + K16161)` (soluble `mmoXYBZ`).
5. **`P_AROM_OXY` — Aerobic Ring-Hydroxylating / Catechol Dioxygenases ($\text{O}_2$-obligate):**
   `(K16243 + K16244 + K16242)` (`dmp/tmo` phenol/toluene monooxygenase) `, K03381` (`catA` catechol 1,2-dioxygenase) `, K00446` (`catE` catechol 2,3-dioxygenase).
6. **`P_THIO_AER` — Aerobic Thiazole Imine Synthesis (`thiO` FAD glycine oxidase, $\text{O}_2$-obligate):**
   `K03153 + K03149` (`thiO` `[EC:1.4.3.19]` producing $\text{H}_2\text{O}_2$ from $\text{O}_2$, paired with `thiG`).

### B. Strictly $\text{O}_2$-Inactivated Metalloenzyme & Radical-SAM Anaerobic Modules ($\mathcal{O}_-$)
Every module in $\mathcal{O}_-$ contains solvent-exposed low-potential `[4Fe-4S]`, `MoFe`/`FeFe`/`Ni-Fe-S` clusters, or Co(I)/radical-SAM intermediates irreversibly inactivated or cleaved by molecular $\text{O}_2$:

1. **`A_NIF` — MoFe Nitrogenase ($\text{O}_2$-inactivated $\text{N}_2 \to \text{NH}_3$):**
   `K02588 + K02586 + K02591` (`nifH + nifD + nifK`, all 3 catalytic subunits required).
2. **`A_WL` — Wood-Ljungdahl CODH/ACS Complex (extreme $\text{O}_2$ sensitivity):**
   `K00198 + (K14138 + K00197 + K00194)` (`cooS/acsA` anaerobic CODH + `acsBCD` acetyl-CoA synthase).
3. **`A_DSR` — Dissimilatory Sulfate Reduction ($\text{SO}_4^{2-} \to \text{H}_2\text{S}$):**
   `(K00394 + K00395) + (K11180 + K11181)` (`aprAB + dsrAB`).
4. **`A_MCR` — Methanogenesis Methyl-CoM Reductase ($\text{O}_2$-inactivated Ni-$F_{430}$):**
   `K00399 + K00401 + K00402` (`mcrA + mcrB + mcrG`).
5. **`A_FEFE_HYD` — Anaerobic `[FeFe]` Hydrogenase ($\text{O}_2$-destroyed H-cluster):**
   `K00533 + K00534` (`hydA + hydB`).
6. **`A_NOS_NOR` — Suboxic/Anoxic Denitrification Terminal Steps ($\text{NO} \to \text{N}_2\text{O} \to \text{N}_2$, $\text{O}_2$-inhibited $\text{Cu}_Z$ center):**
   `(K04561 + K02305)` (`norBC`) and/or `K00376` (`nosZ`).
7. **`A_BZA` — Anaerobic Benzimidazole Lower-Ligand Synthesis (radical-SAM, $\text{O}_2$-inactivated):**
   `(K22465 , K22466)` (`bzaA/B` or `bzaF` 5-hydroxybenzimidazole synthase `[EC:4.1.99.23]`) `+ (K25497 , K25498 , K25499)` (`bzaCDE` methyltransferases).

### C. Biochemically Compatible Vitamin & Amino-Acid Modules (Matched Control Pool $\mathcal{C}_0$)
Redox-neutral cytosolic anabolic pathways compatible across oxic and anoxic cytosols (matched in step count $L = 3\text{–}6$ to the modules above):

1. **`C_HIS` — Histidine biosynthesis ($L=4$ core markers):** `K00765 + K01814 + K01693 + K00013` (`hisG + hisA + hisB + hisD`).
2. **`C_TRP` — Tryptophan biosynthesis ($L=5$):** `(K01657 + K01658) + K00766 + K01609 + (K01695 + K01696)` (`trpEG + trpD + trpC + trpAB`).
3. **`C_LEU` — Leucine biosynthesis ($L=3$):** `K01649 + (K01703 + K01704) + K00052` (`leuA + leuCD + leuB`).
4. **`C_ILV` — Valine/Isoleucine biosynthesis ($L=4$):** `(K01652 , K01653) + K00053 + K01687 + K00826` (`ilvBH + ilvC + ilvD + ilvE`).
5. **`C_LYS_DAP` — Lysine DAP biosynthesis ($L=3$ markers):** `K01714 + K00215 + K01586` (`dapA + dapB + lysA`).
6. **`C_ARG` — Arginine urea-cycle branch ($L=3$, paired with `M00028` ornithine):** `K00611 + K01940 + K01755` (`argF + argG + argH`).
7. **`C_PRO` — Proline biosynthesis ($L=3$):** `K00931 + K00147 + K00286` (`proB + proA + proC`).
8. **`C_RIB` — Riboflavin biosynthesis ($L=4$):** `K01497 + K14652 + K00794 + K00793` (`ribA + ribBA + ribH + ribE`).
9. **`C_PAN` — Pantothenate biosynthesis ($L=3$):** `K00606 + K00077 + K01918` (`panB + panE + panC`).
10. **`C_BIO` — Biotin ring closure ($L=3$):** `K00652 + K00833 + K01012` (`bioF + bioA + bioB`).
11. **`C_THI_SPLIT` — Thiamin two-branch synthesis & precursor salvage:**
    - HMP branch: `K03147 + K00941` (`thiC + thiD`);
    - HET branch: `K03149 + (K03150 , K03153)` (`thiG + (thiH OR thiO)`);
    - Coupling: `K00788 + K00946` (`thiE + thiL`).
12. **`C_MET_B12` — Methionine synthase dependency switch (recipient readout):**
    - Obligate B12-dependent Met user: `K00548` (`metH`) `AND NOT K00549` (`metE` B12-independent synthase absent).

---

## 3. Primary Endpoint Moved to Recipients (Opus `248920882680348672` + Astra `248920851231457280`, `248920953694109696`)

### Why Donor-Only Splitting Cannot Be Primary
Aerobes and anaerobes are distinct organisms across micro-gradients, so $\mathcal{O}_+$ and $\mathcal{O}_-$ pathways sit in disjoint genomes regardless of whether any recipient depends on both (Opus `248920882680348672`, Astra `248920851231457280`). Donor-level partitioning (`co-encoded_only`, `split_only`, `both_present`, `unresolved` per event $\times$ pair, keeping all four in the denominator) is therefore recorded strictly as a **descriptive background check on functional partitioning**, NOT as the test of plexitrophy or blocked helper stacking.

### Primary Recipient-Level Contrast (Falsifiable in Both Directions)
1. **Non-circular recipient eligibility (defined *before* inspecting co-occurring donors, per Astra `248920953694109696`):**
   - Select high-completeness focal MAGs $X$ ($\ge 90\%$ CheckM2, $\le 5\%$ contamination) that lack a specified pair of nutritional endpoints $(N_1, N_2)$ based **solely on $X$'s own genome** (calibrated against the contig-block downsampling null).
2. **Complete single-donor & external-medium alternatives explicitly allowed to win ($d \le 1$, per Astra `248920953694109696`):**
   - For every eligible recipient $X$ in an independent sample/event under specified external supply $E$, enumerate all minimal feasible supports $\mathcal{M}(X, E)$ over co-occurring MAGs (plus `unbinned`).
   - Crucially, **if a single co-occurring donor MAG encodes the complete pathway(s) satisfying $(N_1, N_2)$** — for example, in **I6**, a single bacterium encoding *both* corrinoid synthesis and DMB synthesis (`P_B12_AER + P_BLUB` OR anaerobic `cbiK/X+cbiG + A_BZA`), or in **I2**, dissolved inorganic $\text{NH}_4^+/\text{NO}_3^-$ in $E$ satisfying fixed-N demand without requiring a local `nifHDK` diazotroph — then **$d_S(X) \le 1$ wins**, and $X$ is **NOT** scored as two-donor plexitrophic! Where ambient dissolved metabolite concentrations (e.g., dissolved B12, DIN, DMSP) are unmeasured, results are reported explicitly as **conditional on the stated supply scenario $E$** (testing sensitivity across $E_{\text{minimal}}$ vs. $E_{\text{+DIN}}$ vs. $E_{\text{+ambient\_vitamins}}$).
   - Every eligible recipient $\times$ event is classified into `{d=0 (ambient E suffices), d=1 (single-helper stacked donor suffices), d>=2_or_c>=3 (plexitrophic across all minimal supports), unresolved (empty support, unbinned-only donor, or intermediate step completeness)}`, **retaining `d=0`, `d=1`, and `unresolved` in the reported denominator**.
3. **Matched Contrast (Lead pairs I6, I3, plus split cofactor pathways vs. compatible pairs in $\mathcal{C}_0 \times \mathcal{C}_0$):**
   - Compare the probability $P(\text{plexitrophic: } \forall S,\, d_S \ge 2 \lor c_S \ge 3 \mid X \text{ eligible})$ between recipients missing an **incompatible-producer pair** (lead pair **I6**: `metH`-only `K00548 \ K00549` lacking corrinoid ring AND lower-ligand synthesis, where complete B12 requires corrinoid ring + DMB lower ligand; **I3**: lacking B12 AND assimilatory sulfate reduction `M00176` where reduced S + B12 span redox zones; **I2**: under explicit DIN-depleted $E_{\text{minimal}}$ only) versus recipients missing a **matched compatible pair** from $\mathcal{C}_0 \times \mathcal{C}_0$ (e.g., missing two redox-neutral amino acids/vitamins such as `C_HIS + C_TRP` or `C_LEU + C_PAN`), matched on `(total_missing_modules_in_X, genome_size_X, CheckM2_completeness_X)` with crossed random intercepts for event, function identity, and GTDB clade.
   - Per Opus (`248920882680348672` point 3), **I1 (`P_PSII × A_NIF`)** is dominated by compartmentalizing/circadian cyanobacteria and is reported **separately** rather than pooled into the recipient estimate.

| Pair ID | Function 1 ($\mathcal{O}_+$: $\text{O}_2$-evolving / $\text{O}_2$-obligate) | Function 2 ($\mathcal{O}_-$: $\text{O}_2$-inactivated / anaerobic) | Recipient $X$ Eligibility (Genome of $X$ Only) & $d=1$ Single-Donor Falsifier | Known Single-Cell Bypass to Tally Explicitly |
|---|---|---|---|---|
| **I1** | `P_PSII` (oxygenic photosynthesis) | `A_NIF` (`nifHDK` $\text{N}_2$ fixation) | Fixed C + fixed N for C/N-auxotrophic heterotrophs | Cyanobacterial heterocysts (`K25172` `hgd`) or circadian temporal separation (`K08480+K08481+K08482` `kaiABC` in unicellular diazotrophs) |
| **I2** | `P_BLUB` (`bluB` $\text{O}_2$-obligate DMB lower ligand) | `A_NIF` (`nifHDK` $\text{N}_2$ fixation) or `A_WL` (`CODH/ACS`) | Cobalamin (`metH`-only `K00548 \ K00549` recipients) + fixed N / acetate across oxic–anoxic microniches | Respiratory $\text{O}_2$ scavenging via high-affinity cytochrome *bd* (`K00425+K00426` `cydAB`) or *cbb3* (`K00404+K00405` `ccoNO`) |
| **I3** | `P_B12_AER` (`cobG+cobNST` aerobic corrinoid synthesis) | `A_DSR` (`aprAB+dsrAB` $\text{H}_2\text{S}$/reduced-S supply) | B12 + reduced sulfur (for assimilatory sulfate-reduction auxotrophs lacking `M00176`, e.g. SAR11-like) | None known within one non-compartmentalized cytosol |
| **I4** | `P_CU_MMO` (`amoABC`/`pmoABC` aerobic N/CH4 oxidation) | `A_NOS_NOR` (`norBC`/`nosZ` suboxic denitrification) | Coupled nitrification–denitrification or methane-driven N-loss across particle/biofilm redox gradients | Microaerophilic `cbb3`/`cydAB` respiratory protection ("nitrifier denitrification" carries `nirK`/`norBC`, rarely `nosZ`) |
| **I5** | `P_AROM_OXY` (aerobic aromatic mono/dioxygenases) | `A_DSR` or `A_NOS_NOR` or `A_NIF` | Hydrocarbon/DOM oxidation coupled to N/S supply in oil-snow / sediment / ice-brine aggregates | Facultative anaerobes switching regulons temporally |
| **I6** | `P_BLUB` / `P_B12_AER` (aerobic DMB / B12) | Anaerobic corrinoid ring synthesis (`cbiK/X` `K02190,K03795` + `cbiG` `K13541` without `bluB`) | Complete cobalamin via **inter-species cobinamide + DMB lower-ligand remodeling** (`cobT+cobS+cobC` `K00768+K02233+K02226`) | Single genome encoding anaerobic `bza` (`A_BZA`) alongside anaerobic `cbi` |

---

## 4. Pre-Registered Intracellular Compartmentalization & Respiratory-Protection Accounting (Opus Challenge #1)

To ensure P2 **cannot confirm itself** by ignoring single-cell workarounds (`248917417287663616`), every MAG that co-encodes an incompatible pair $(F_1 \in \mathcal{O}_+,\; F_2 \in \mathcal{O}_-)$ is classified into one of three mutually exclusive buckets (with all counts and denominators reported):

1. **Spatial compartmentalizer ($S_{\text{spatial}}$):** Co-encodes heterocyst glycolipid marker (`K25172`) or bacterial microcompartment shell proteins (`K27260`, `K27263`, `K27266`–`K27272`, `K04027` `eutM`).
2. **Temporal / respiratory-protection specialist ($S_{\text{temp/resp}}$):** Co-encodes circadian oscillator (`K08480 + K08481 + K08482` `kaiABC` in non-heterocystous cyanobacteria) or high-affinity $\text{O}_2$-scavenging terminal oxidases (`K00425 + K00426` `cydAB` cytochrome *bd* or `K00404 + K00405` `ccoNO` *cbb3*-type oxidase, classic *Azotobacter* respiratory protection).
3. **Unexplained single-genome stacker ($S_{\text{unexplained}}$):** Co-encodes $(F_1, F_2)$ on contamination-screened ($\le 5\%$ CheckM2 contamination, concordant contig coverage/tetranucleotide profile) contigs **without** any marker in (1) or (2). *Note:* any apparent co-encoder with discordant contig depth or $>5\%$ CheckM2 contamination is flagged as a **chimeric-bin candidate** rather than counted as a biological stacker.

### Decision Bands Around a Single Frozen Margin $\delta_{\text{OR}}$ (Recipient-Level Contrast, Zero Gaps)
Let $[L, U]$ be the 95% confidence/credible interval (cluster-bootstrapped over independent events and function identities, matched on recipient missing-module count, genome size, CheckM2 completeness, and GTDB clade) for the recipient-level log-odds difference
$$\Delta \theta_{\text{recip}}(E) = \text{logit}\, P\big(d \ge 2 \lor c \ge 3 \mid X \in \text{Incompatible-pair auxotrophs},\, E\big) - \text{logit}\, P\big(d \ge 2 \lor c \ge 3 \mid X \in \mathcal{C}_0 \times \mathcal{C}_0 \text{ auxotrophs},\, E\big),$$
evaluated against a single pre-registered meaningful-effect margin $\delta_{\text{OR}}$ (e.g., $\delta_{\text{OR}} = \ln 2 \approx 0.69$, fixed in `protocols/` before inspecting data):
- **Supported (blocked helper stacking increases recipient plexitrophy under $E$):** $L > \delta_{\text{OR}}$ (entire 95% CI strictly above $\delta_{\text{OR}}$).
- **Falsified (single-helper stacking or complete-pathway donors satisfy incompatible-pair recipients as often as compatible-pair recipients):** $U < \delta_{\text{OR}}$ (entire 95% CI strictly below $\delta_{\text{OR}}$).
- **Inconclusive:** $L \le \delta_{\text{OR}} \le U$ (95% CI straddles $\delta_{\text{OR}}$), OR $< 10$ independent events contain eligible high-completeness recipients in both strata.

Donor-level partitioning (`co-encoded_only`, `split_only`, `both_present`, `unresolved`) and single-genome co-encoder bypass categories $(S_{\text{spatial}}, S_{\text{temp/resp}}, S_{\text{unexplained}})$ are reported alongside as descriptive mechanism accounting, without creating a second conflicting decision rule.
