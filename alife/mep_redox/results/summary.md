## T1. Entropy production (J K^-1 m^-2 d^-1, mean ± sd over seeds; last 300 d)

| scenario | n | sigma_chemotroph | NPP (mmol C m^-2 d^-1) | sigma_chem per NPP (J K^-1 per mmol C) | sigma_algae | algal share | fixed-enzyme copy (gauge) | greedy (H->0) | random alloc (mean) | best-found ref (steady) | emergent / best-found | gauge rel. diff | fixed-field sigma within 0.5% | Powell converged | 2nd-pass rel. gain | params at bound |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| neutral | 3 | 56.2 ± 26.5 | 54.5 ± 26.4 | 1.039 ± 0.021 | 219.5 ± 106.2 | 0.795 ± 0.003 | 56.3 ± 26.5 | 55.7 ± 26.1 | 56.3 ± 26.7 | 58.4 ± 27.5 | 0.962 ± 0.002 | 0.0003 ± 0.0004 | 3/3 | 0/3 | 0.0003 ± 0.0001 | 0.0 ± 0.0 |
| crowd | 3 | 45.3 ± 19.1 | 43.4 ± 18.9 | 1.047 ± 0.014 | 174.9 ± 76.1 | 0.794 ± 0.002 | 45.3 ± 19.1 | 44.8 ± 18.7 | 44.9 ± 19.2 | 47.0 ± 19.8 | 0.963 ± 0.003 | 0.0001 ± 0.0002 | 3/3 | 0/3 | 0.0005 ± 0.0001 | 0.0 ± 0.0 |
| stackadv | 3 | 58.4 ± 28.1 | 56.6 ± 27.9 | 1.039 ± 0.021 | 228.1 ± 112.5 | 0.795 ± 0.003 | 58.4 ± 28.1 | 57.9 ± 27.7 | 58.7 ± 28.3 | 60.6 ± 29.1 | 0.962 ± 0.002 | 0.0003 ± 0.0002 | 3/3 | 0/3 | 0.0003 ± 0.0001 | 0.0 ± 0.0 |
| constitutive | 3 | 109.4 ± 21.6 | 108.9 ± 22.2 | 1.006 ± 0.007 | 438.0 ± 89.4 | 0.800 ± 0.001 | 109.5 ± 21.6 | 107.5 ± 20.9 | 107.3 ± 21.0 | 113.5 ± 22.4 | 0.964 ± 0.002 | 0.0008 ± 0.0003 | 3/3 | 0/3 | 0.0006 ± 0.0001 | 0.0 ± 0.0 |
| channel | 3 | 56.3 ± 26.5 | 54.4 ± 26.4 | 1.041 ± 0.020 | 219.4 ± 106.1 | 0.795 ± 0.003 | 56.2 ± 26.4 | 55.7 ± 26.0 | 56.2 ± 26.6 | 58.4 ± 27.4 | 0.964 ± 0.003 | 0.0021 ± 0.0012 | 3/3 | 0/3 | 0.0003 ± 0.0001 | 0.0 ± 0.0 |
| toxic_channel | 3 | 52.3 ± 23.8 | 50.4 ± 23.5 | 1.043 ± 0.018 | 203.3 ± 94.6 | 0.794 ± 0.003 | 52.2 ± 23.7 | 51.8 ± 23.3 | 52.3 ± 23.9 | 54.2 ± 24.5 | 0.964 ± 0.003 | 0.0021 ± 0.0012 | 3/3 | 0/3 | 0.0004 ± 0.0001 | 0.0 ± 0.0 |
| electron_yield | 3 | 39.9 ± 17.5 | 38.0 ± 16.9 | 1.055 ± 0.011 | 152.8 ± 68.2 | 0.792 ± 0.002 | 39.9 ± 17.5 | 39.2 ± 17.1 | 38.6 ± 17.6 | 41.3 ± 18.2 | 0.968 ± 0.002 | 0.0002 ± 0.0001 | 3/3 | 0/3 | 0.0009 ± 0.0005 | 0.0 ± 0.0 |
| dispersal | 3 | 49.4 ± 21.8 | 47.6 ± 21.7 | 1.043 ± 0.018 | 191.4 ± 87.2 | 0.794 ± 0.003 | 49.4 ± 21.9 | 48.5 ± 21.2 | 48.2 ± 21.2 | 51.3 ± 22.6 | 0.962 ± 0.002 | 0.0007 ± 0.0006 | 3/3 | 0/3 | 0.0001 ± 0.0002 | 0.0 ± 0.0 |
| forced | 3 | 45.3 ± 23.3 | 43.1 ± 22.6 | 1.055 ± 0.015 | 175.2 ± 91.8 | 0.794 ± 0.002 | na | na | na | na | na | na | 0/0 | 0/0 | na | na |

## T2. Packaging contrast: same system, pool restricted to split or stacked genomes

| scenario | sigma full pool | sigma split-only | sigma stack-only | split vs stack (%) | per-NPP full | per-NPP split-only | per-NPP stack-only | per-NPP split vs stack (%) |
|---|---|---|---|---|---|---|---|---|
| neutral | 56.2 ± 26.5 | 88.0 ± 54.2 | 71.8 ± 35.9 | 16.14 ± 21.45 | 1.039 ± 0.021 | 1.022 ± 0.032 | 1.030 ± 0.027 | -0.80 ± 0.46 |
| crowd | 45.3 ± 19.1 | 54.9 ± 32.1 | 52.5 ± 24.3 | 1.00 ± 13.44 | 1.047 ± 0.014 | 1.040 ± 0.024 | 1.039 ± 0.021 | 0.02 ± 0.42 |
| stackadv | 58.4 ± 28.1 | 106.4 ± 65.5 | 75.9 ± 39.3 | 32.88 ± 25.20 | 1.039 ± 0.021 | 1.015 ± 0.035 | 1.030 ± 0.029 | -1.46 ± 0.57 |
| constitutive | 109.4 ± 21.6 | 109.3 ± 21.7 | 111.4 ± 19.3 | -2.25 ± 3.50 | 1.006 ± 0.007 | 1.005 ± 0.007 | 0.998 ± 0.007 | 0.74 ± 0.19 |
| channel | 56.3 ± 26.5 | 88.0 ± 54.2 | 72.0 ± 36.1 | 15.85 ± 21.27 | 1.041 ± 0.020 | 1.022 ± 0.032 | 1.032 ± 0.027 | -0.95 ± 0.54 |
| toxic_channel | 52.3 ± 23.8 | 80.1 ± 49.0 | 66.0 ± 32.3 | 15.14 ± 22.40 | 1.043 ± 0.018 | 1.023 ± 0.031 | 1.034 ± 0.025 | -1.07 ± 0.65 |

## T3. Division of labour (full pool). stacked share = fraction of consumer-step flux done by genomes that also carry the producer step; depth_sep = consumer minus producer flux-weighted depth (m)

| scenario | stacked: denitrification NO3->NO2->N2 | stacked: nitrification NH4->NO2->NO3 | stacked: fermentation->methanogenesis (H2) | sep: denitrification NO3->NO2->N2 | sep: nitrification NH4->NO2->NO3 | sep: fermentation->methanogenesis (H2) | sep: fermentation->H2 sulfate reduction | sep: sulfate reduction->sulfide oxidation | sep: methanogenesis->methanotrophy |
|---|---|---|---|---|---|---|---|---|---|
| neutral | 0.99 ± 0.00 | 0.81 ± 0.06 | 0.14 ± 0.19 | 0.5 ± 0.1 | 1.6 ± 0.5 | -0.0 ± 0.1 | -0.1 ± 0.1 | -12.5 ± 2.1 | -9.0 ± 2.7 |
| crowd | 0.93 ± 0.05 | 0.12 ± 0.07 | 0.13 ± 0.20 | 0.4 ± 0.1 | 1.7 ± 0.4 | -0.0 ± 0.1 | -0.2 ± 0.2 | -12.9 ± 3.2 | -8.4 ± 2.1 |
| stackadv | 1.00 ± 0.00 | 0.94 ± 0.03 | 0.16 ± 0.18 | 0.4 ± 0.1 | 1.5 ± 0.6 | -0.0 ± 0.1 | -0.1 ± 0.1 | -12.7 ± 2.5 | -9.2 ± 2.8 |
| constitutive | 0.25 ± 0.16 | 0.92 ± 0.08 | 0.42 ± 0.28 | 0.5 ± 0.3 | 0.8 ± 0.4 | 0.6 ± 0.5 | -0.3 ± 0.2 | -11.3 ± 0.3 | -12.6 ± 1.1 |
| channel | 1.00 ± 0.00 | 0.96 ± 0.02 | 0.26 ± 0.29 | 0.5 ± 0.1 | 1.4 ± 0.8 | -0.0 ± 0.1 | -0.1 ± 0.1 | -12.5 ± 2.1 | -9.0 ± 2.7 |
| toxic_channel | 0.99 ± 0.00 | 0.96 ± 0.03 | 0.25 ± 0.29 | 0.3 ± 0.1 | 1.4 ± 0.7 | -0.0 ± 0.1 | -0.1 ± 0.1 | -13.5 ± 3.7 | -8.9 ± 2.5 |
| electron_yield | 0.96 ± 0.03 | 0.28 ± 0.14 | 0.03 ± 0.05 | 0.6 ± 0.5 | 4.3 ± 4.1 | 0.1 ± 0.1 | -0.0 ± 0.1 | -15.1 ± 4.8 | -7.3 ± 1.4 |
| dispersal | 1.00 ± 0.00 | 0.87 ± 0.01 | 0.15 ± 0.21 | 0.3 ± 0.1 | 1.8 ± 0.4 | 0.1 ± 0.2 | -0.3 ± 0.1 | -13.8 ± 3.3 | -10.2 ± 2.2 |
| forced | 0.72 ± 0.48 | 0.81 ± 0.05 | 0.13 ± 0.19 | 0.0 ± 0.4 | 2.2 ± 1.6 | -0.0 ± 0.2 | -0.0 ± 0.1 | -14.1 ± 5.5 | -8.2 ± 2.7 |

## T4. Realised partners vs counterfactual necessity (knock-outs, same medium, endpoint = recipient column biomass >= 10% of an UNPERTURBED CONTROL at the same horizon; d_lb = min hitting set of failing removal sets (direct: own-input groups + direct donors; upstream: + indirect), lower bounds conditional on monotone rescue)

| run | recipient | d_real | c_real | n single removals failing | failing: direct donors | failing: indirect | solver failures | d_lb_direct | d_lb_upstream | failing OR-groups | algal share of intake |
|---|---|---|---|---|---|---|---|---|---|---|---|
| neutral__s0 | AER | 1 | 9 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG | 1.00 |
| neutral__s0 | ROS | 1 | 9 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG;DMSP:ALG | 1.00 |
| neutral__s0 | AERNAR | 3 | 9 | 2 | ALG | ROS | 0 | 1 | 2 | DOC:ALG;O2:ALG | 0.99 |
| neutral__s0 | AERDEN | 4 | 9 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG | 0.96 |
| neutral__s0 | SRBX | 3 | 9 | 2 | ALG | ROS | 0 | 1 | 2 | DOC:ALG | 0.94 |
| neutral__s0 | AOB | 6 | 9 | 2 | ROS,ALG |  | 0 | 2 | 2 | O2:ALG;NH4:AER|AERDEN|AERNAR|ROS|SRBX | 0.18 |
| neutral__s0 | CMX | 7 | 9 | 1 | ALG |  | 0 | 2 | 2 | O2:ALG;NH4:AER|AERDEN|AERNAR|ROS|SRBX | 0.14 |
| neutral__s0 | AMX | 5 | 1 | 2 | ROS | ALG | 0 | 1 | 2 | NH4:AER|AERDEN|AERNAR|ROS|SRBX | 0.00 |
| neutral__s0 | SOB | 3 | 1 | 2 | ROS,ALG |  | 0 | 2 | 2 | O2:ALG;H2S:ROS|SRBX | 0.32 |
| neutral__s0 | SOBX | 5 | 9 | 2 | ROS,ALG |  | 0 | 2 | 2 | O2:ALG;H2S:ROS|SRBX | 0.20 |
| neutral__s0 | ALG | 4 | 9 | 0 |  |  | 0 | 0 | 0 |  | 0.00 |
| neutral__s1 | AER | 1 | 8 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG | 1.00 |
| neutral__s1 | ROS | 1 | 8 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG;DMSP:ALG | 1.00 |
| neutral__s1 | AERNAR | 3 | 8 | 2 | ALG | ROS | 0 | 1 | 2 | DOC:ALG;O2:ALG | 0.97 |
| neutral__s1 | AERDEN | 5 | 8 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG | 0.95 |
| neutral__s1 | SRBX | 3 | 8 | 2 | ALG | ROS | 0 | 1 | 2 | DOC:ALG | 0.96 |
| neutral__s1 | AOB | 5 | 1 | 1 | ALG |  | 0 | 2 | 2 | O2:ALG;NH4:AER|AERDEN|AERNAR|ROS|SRBX | 0.20 |
| neutral__s1 | CMX | 5 | 8 | 1 | ALG |  | 0 | 2 | 2 | O2:ALG;NH4:AER|AERDEN|AERNAR|ROS|SRBX | 0.15 |
| neutral__s1 | AMX | 4 | 1 | 2 | ROS | ALG | 0 | 1 | 2 | NH4:AER|AERDEN|AERNAR|ROS|SRBX | 0.00 |
| neutral__s1 | SOB | 3 | 1 | 2 | ROS,ALG |  | 0 | 2 | 2 | O2:ALG;H2S:ROS|SRBX | 0.34 |
| neutral__s1 | SOBX | 5 | 8 | 2 | ROS,ALG |  | 0 | 2 | 2 | O2:ALG;H2S:ROS|SRBX | 0.16 |
| neutral__s1 | ALG | 5 | 8 | 0 |  |  | 0 | 0 | 0 |  | 0.00 |
| neutral__s2 | AER | 1 | 6 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG | 1.00 |
| neutral__s2 | ROS | 1 | 6 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG;DMSP:ALG | 1.00 |
| neutral__s2 | AERNAR | 3 | 6 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG | 0.97 |
| neutral__s2 | AERDEN | 4 | 6 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG | 0.95 |
| neutral__s2 | AOB | 4 | 1 | 0 |  |  | 0 | 0 | 0 |  | 0.14 |
| neutral__s2 | CMX | 4 | 6 | 1 | ALG |  | 0 | 2 | 2 | O2:ALG;NH4:AER|AERDEN|AERNAR|ROS | 0.10 |
| neutral__s2 | SOBX | 5 | 1 | 2 | ROS,ALG |  | 0 | 2 | 2 | O2:ALG;H2S:ROS | 0.13 |
| neutral__s2 | ALG | 4 | 6 | 0 |  |  | 0 | 0 | 0 |  | 0.00 |
| crowd__s0 | AER | 1 | 10 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG | 1.00 |
| crowd__s0 | ROS | 1 | 10 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG;DMSP:ALG | 1.00 |
| crowd__s0 | AERNAR | 2 | 10 | 2 | ALG | ROS | 0 | 1 | 2 | DOC:ALG;O2:ALG | 0.95 |
| crowd__s0 | AERDEN | 5 | 10 | 2 | ALG | ROS | 0 | 1 | 2 | DOC:ALG;O2:ALG | 0.94 |
| crowd__s0 | SRB | 1 | 10 | 2 | ALG | ROS | 0 | 1 | 2 | DOC:ALG | 0.97 |
| crowd__s0 | SRBX | 3 | 10 | 2 | ALG | ROS | 0 | 1 | 2 | DOC:ALG | 0.92 |
| crowd__s0 | AOB | 5 | 10 | 1 | ALG |  | 0 | 2 | 2 | O2:ALG;NH4:AER|AERDEN|AERNAR|ROS|SRB|SRBX | 0.20 |
| crowd__s0 | NOB | 2 | 10 | 1 | ALG |  | 0 | 1 | 1 | O2:ALG | 0.10 |
| crowd__s0 | AMX | 5 | 1 | 2 | ROS | ALG | 0 | 1 | 2 | NH4:AER|AERDEN|AERNAR|ROS|SRB|SRBX | 0.00 |
| crowd__s0 | SOB | 4 | 1 | 2 | ROS,ALG |  | 0 | 2 | 2 | O2:ALG;H2S:ROS|SRB|SRBX | 0.34 |
| crowd__s0 | SOBX | 4 | 10 | 2 | ROS,ALG |  | 0 | 2 | 2 | O2:ALG;H2S:ROS|SRB|SRBX | 0.07 |
| crowd__s0 | ALG | 5 | 10 | 0 |  |  | 0 | 0 | 0 |  | 0.00 |
| constitutive__s0 | AER | 1 | 3 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG | 1.00 |
| constitutive__s0 | ROS | 1 | 3 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG;DMSP:ALG | 1.00 |
| constitutive__s0 | AERNAR | 3 | 1 | 2 | ALG | ROS | 0 | 1 | 2 | DOC:ALG;O2:ALG | 0.88 |
| constitutive__s0 | SRB | 2 | 2 | 2 | ALG | ROS | 0 | 1 | 2 | DOC:ALG | 0.95 |
| constitutive__s0 | SRBX | 3 | 1 | 2 | ALG | ROS | 0 | 1 | 2 | DOC:ALG | 0.31 |
| constitutive__s0 | AOB | 2 | 1 | 2 | ROS,ALG |  | 0 | 2 | 2 | O2:ALG;NH4:AER|ROS|SRB | 0.11 |
| constitutive__s0 | CMX | 3 | 1 | 2 | ROS,ALG |  | 0 | 2 | 2 | O2:ALG;NH4:AER|ROS|SRB | 0.06 |
| constitutive__s0 | AMX | 3 | 1 | 1 |  | ALG | 0 | 0 | 1 |  | 0.00 |
| constitutive__s0 | SOB | 4 | 2 | 1 | ALG |  | 0 | 1 | 1 | O2:ALG | 0.21 |
| constitutive__s0 | SOBX | 5 | 1 | 2 | ROS,ALG |  | 0 | 2 | 2 | O2:ALG;H2S:ROS|SRB|SRBX | 0.05 |
| constitutive__s0 | MOB | 3 | 1 | 0 |  |  | 0 | 0 | 0 |  | 0.23 |
| constitutive__s0 | ALG | 2 | 3 | 0 |  |  | 0 | 0 | 0 |  | 0.00 |
| channel__s0 | AER | 1 | 8 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG | 1.00 |
| channel__s0 | ROS | 1 | 8 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG;DMSP:ALG | 1.00 |
| channel__s0 | AERNAR | 2 | 8 | 2 | ALG | ROS | 0 | 1 | 2 | DOC:ALG;O2:ALG | 0.98 |
| channel__s0 | AERDEN | 3 | 8 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG | 0.96 |
| channel__s0 | SRBX | 3 | 8 | 2 | ALG | ROS | 0 | 1 | 2 | DOC:ALG | 0.94 |
| channel__s0 | CMX | 6 | 8 | 1 | ALG |  | 0 | 2 | 2 | O2:ALG;NH4:AER|AERDEN|AERNAR|ROS|SRBX | 0.18 |
| channel__s0 | AMX | 4 | 1 | 2 | ROS | ALG | 0 | 1 | 2 | NH4:AER|AERDEN|AERNAR|ROS|SRBX | 0.00 |
| channel__s0 | SOB | 3 | 1 | 2 | ROS,ALG |  | 0 | 2 | 2 | O2:ALG;H2S:ROS|SRBX | 0.32 |
| channel__s0 | SOBX | 4 | 8 | 2 | ROS,ALG |  | 0 | 2 | 2 | O2:ALG;H2S:ROS|SRBX | 0.20 |
| channel__s0 | ALG | 4 | 8 | 0 |  |  | 0 | 0 | 0 |  | 0.00 |
| electron_yield__s0 | AER | 1 | 11 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG | 1.00 |
| electron_yield__s0 | ROS | 1 | 11 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG;DMSP:ALG | 1.00 |
| electron_yield__s0 | AERNAR | 3 | 1 | 2 | ALG | ROS | 0 | 1 | 2 | DOC:ALG;O2:ALG | 0.99 |
| electron_yield__s0 | AERDEN | 4 | 11 | 1 | ALG |  | 0 | 1 | 1 | DOC:ALG;O2:ALG | 0.98 |
| electron_yield__s0 | SRBH | 2 | 11 | 2 |  | ROS,ALG | 0 | 1 | 3 | H2:FER|FERMET | 0.00 |
| electron_yield__s0 | SRBX | 3 | 1 | 0 |  |  | 0 | 0 | 0 |  | 0.64 |
| electron_yield__s0 | FER | 1 | 11 | 2 | ALG | ROS | 0 | 1 | 2 | DOC:ALG | 1.00 |
| electron_yield__s0 | MET | 2 | 1 | 2 |  | ROS,ALG | 0 | 1 | 3 | H2:FER|FERMET | 0.00 |
| electron_yield__s0 | FERMET | 2 | 11 | 2 | ALG | ROS | 0 | 1 | 2 | DOC:ALG | 0.91 |
| electron_yield__s0 | AOB | 4 | 11 | 1 | ALG |  | 0 | 2 | 2 | O2:ALG;NH4:AER|AERDEN|AERNAR|FER|FERMET|ROS | 0.20 |
| electron_yield__s0 | NOB | 4 | 11 | 1 | ALG |  | 0 | 2 | 2 | O2:ALG;NO2:AOB|CMX|SOBX | 0.09 |
| electron_yield__s0 | CMX | 6 | 11 | 1 | ALG |  | 0 | 2 | 2 | O2:ALG;NH4:AER|AERDEN|AERNAR|FER|FERMET|ROS | 0.19 |
| electron_yield__s0 | AMX | 3 | 1 | 2 | ROS | ALG | 0 | 1 | 2 |  | 0.00 |
| electron_yield__s0 | SOB | 3 | 1 | 2 | ROS,ALG |  | 0 | 2 | 2 | O2:ALG;H2S:ROS|SRBH|SRBX | 0.35 |
| electron_yield__s0 | SOBX | 5 | 11 | 2 | ROS,ALG |  | 0 | 2 | 2 | O2:ALG;H2S:ROS|SRBH|SRBX | 0.19 |
| electron_yield__s0 | MOB | 2 | 1 | 2 | ALG | ROS | 0 | 2 | 3 | O2:ALG;CH4:FERMET|MET | 0.38 |
| electron_yield__s0 | ALG | 4 | 11 | 0 |  |  | 0 | 0 | 0 |  | 0.00 |

### T4b. Realised d>=2 vs required d>=2 (recipients with >=0.2% biomass)

| run | recipients | realised d_real>=2 | d_lb_direct>=2 | d_lb_upstream>=2 | realised>=2 but d_lb_direct<2 |
|---|---|---|---|---|---|
| neutral__full__s0 | 11 | 9 | 4 | 7 | 5 |
| neutral__full__s1 | 11 | 9 | 4 | 7 | 5 |
| neutral__full__s2 | 8 | 6 | 2 | 2 | 4 |
| crowd__full__s0 | 12 | 9 | 3 | 8 | 6 |
| constitutive__full__s0 | 12 | 10 | 3 | 6 | 7 |
| channel__full__s0 | 10 | 8 | 3 | 6 | 5 |
| electron_yield__full__s0 | 17 | 14 | 6 | 12 | 8 |

## T5. Observer-side signatures and inference on synthetic data (by configuration)

| config | n | sigma_chem | rho RNA~power | rho RNA~Monod | rho RNA~DNA | unexploited power | KO-rule precision | KO-rule recall | redox path (inferred) | redox path (true) | co-occ precision | co-occ depth-adjacent false edges |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| emergent:neutral | 3 | 56.2 ± 26.5 | 0.76 ± 0.06 | 0.68 ± 0.05 | 0.52 ± 0.07 | 0.053 ± 0.053 | 0.09 ± 0.04 | 0.92 ± 0.13 | 1.45 ± 0.02 | 1.45 ± 0.39 | 0.08 ± 0.02 | 15.0 ± 0.0 |
| best_found_ref | 24 | 60.6 ± 29.4 | 0.69 ± 0.07 | 0.56 ± 0.11 | 0.66 ± 0.11 | 0.056 ± 0.047 | 0.33 ± 0.06 | 0.71 ± 0.24 | 1.36 ± 0.17 | 1.31 ± 0.16 | 0.29 ± 0.05 | 4.8 ± 1.4 |
| greedy | 24 | 57.6 ± 27.8 | 0.84 ± 0.04 | 0.70 ± 0.05 | 0.97 ± 0.01 | 0.140 ± 0.038 | 0.00 ± 0.00 | na | 1.23 ± 0.12 | na | 0.00 ± 0.00 | 6.1 ± 1.1 |
| random1 | 24 | 57.5 ± 27.5 | 0.73 ± 0.05 | 0.56 ± 0.06 | 0.47 ± 0.15 | 0.023 ± 0.024 | 0.45 ± 0.03 | 0.79 ± 0.15 | 1.43 ± 0.10 | 1.52 ± 0.03 | na | 0.0 ± 0.0 |
| random2 | 24 | 58.6 ± 28.6 | 0.73 ± 0.05 | 0.51 ± 0.12 | 0.43 ± 0.17 | 0.015 ± 0.020 | 0.43 ± 0.03 | 0.80 ± 0.15 | 1.43 ± 0.10 | 1.48 ± 0.13 | na | 0.0 ± 0.0 |
| emergent:crowd | 3 | 45.3 ± 19.1 | 0.77 ± 0.08 | 0.67 ± 0.05 | 0.58 ± 0.09 | 0.046 ± 0.077 | 0.11 ± 0.05 | 0.82 ± 0.11 | 1.45 ± 0.02 | 1.39 ± 0.34 | 0.13 ± 0.06 | 14.3 ± 1.5 |
| emergent:stackadv | 3 | 58.4 ± 28.1 | 0.76 ± 0.04 | 0.67 ± 0.05 | 0.51 ± 0.09 | 0.066 ± 0.053 | 0.07 ± 0.03 | 0.91 ± 0.16 | 1.45 ± 0.02 | 1.48 ± 0.42 | 0.05 ± 0.02 | 15.3 ± 1.2 |
| emergent:constitutive | 3 | 109.4 ± 21.6 | 0.80 ± 0.00 | 0.70 ± 0.01 | 0.57 ± 0.04 | 0.026 ± 0.032 | 0.08 ± 0.04 | 0.84 ± 0.14 | 1.45 ± 0.02 | 1.25 ± 0.22 | 0.08 ± 0.04 | 14.3 ± 1.5 |
| emergent:channel | 3 | 56.3 ± 26.5 | 0.76 ± 0.06 | 0.69 ± 0.06 | 0.51 ± 0.09 | 0.048 ± 0.059 | 0.08 ± 0.03 | 0.92 ± 0.13 | 1.45 ± 0.02 | 1.46 ± 0.40 | 0.06 ± 0.03 | 15.7 ± 1.5 |
| emergent:toxic_channel | 3 | 52.3 ± 23.8 | 0.76 ± 0.08 | 0.69 ± 0.03 | 0.53 ± 0.06 | 0.047 ± 0.059 | 0.07 ± 0.03 | 0.89 ± 0.12 | 1.45 ± 0.02 | 1.42 ± 0.36 | 0.04 ± 0.02 | 15.7 ± 0.6 |
| emergent:electron_yield | 3 | 39.9 ± 17.5 | 0.69 ± 0.09 | 0.65 ± 0.08 | 0.64 ± 0.05 | 0.047 ± 0.077 | 0.22 ± 0.17 | 0.94 ± 0.10 | 1.39 ± 0.08 | 1.48 ± 0.29 | 0.25 ± 0.19 | 12.3 ± 5.1 |
| emergent:dispersal | 3 | 49.4 ± 21.8 | 0.73 ± 0.06 | 0.65 ± 0.03 | 0.44 ± 0.13 | 0.193 ± 0.068 | 0.19 ± 0.07 | 0.74 ± 0.25 | 1.00 ± 0.00 | 1.39 ± 0.34 | 0.26 ± 0.06 | 8.0 ± 1.7 |
| emergent:forced | 3 | 45.3 ± 23.3 | 0.74 ± 0.07 | 0.63 ± 0.06 | 0.54 ± 0.07 | 0.044 ± 0.066 | 0.09 ± 0.03 | 0.90 ± 0.12 | 1.46 ± 0.04 | 1.45 ± 0.39 | 0.08 ± 0.01 | 15.0 ± 1.7 |

## T6. Horizon test (neutral, seed 0; sterile start). Mean chemotrophic sigma over [0, H]

| policy | H=1 d | H=10 d | H=100 d | H=1000 d |
|---|---|---|---|---|
| static_best_H1 | 1877.6 | 999.6 | 180.7 | 78.2 |
| static_best_H10 | 1645.0 | 1092.1 | 207.9 | 91.9 |
| static_best_H100 | 984.9 | 763.3 | 218.5 | 94.4 |
| static_best_H1000 | 917.6 | 684.9 | 216.5 | 94.8 |
| emergent_field_static | 1519.7 | 1083.6 | 211.7 | 92.9 |
| greedy_instantaneous | 1870.7 | 1200.3 | 214.8 | 91.9 |
| population_assembly | 11.2 | 1136.1 | 252.7 | 97.2 |

## T7. Column biomass by type (mmol C m^-2, mean over seeds)

| scenario | AER | ROS | NAR | NIR | DEN | AERNAR | AERDEN | SRB | SRBH | SRBX | FER | MET | FERMET | AOB | NOB | CMX | AMX | SOB | SNR | SOBX | MOB | ALG |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| neutral | 82 | 62 | 0 | 0 | 0 | 65 | 95 | 0 | 0 | 8 | 0 | 0 | 0 | 4 | 1 | 21 | 2 | 24 | 0 | 11 | 1 | 231 |
| crowd | 127 | 48 | 0 | 0 | 0 | 42 | 36 | 3 | 0 | 3 | 0 | 0 | 0 | 16 | 3 | 2 | 2 | 20 | 1 | 2 | 1 | 198 |
| stackadv | 50 | 65 | 0 | 0 | 0 | 66 | 124 | 0 | 0 | 9 | 0 | 0 | 0 | 1 | 0 | 27 | 2 | 14 | 0 | 24 | 1 | 237 |
| constitutive | 237 | 123 | 0 | 0 | 0 | 21 | 0 | 22 | 0 | 2 | 1 | 0 | 0 | 27 | 1 | 18 | 9 | 74 | 2 | 6 | 2 | 374 |
| channel | 81 | 62 | 0 | 0 | 0 | 65 | 96 | 0 | 0 | 8 | 0 | 0 | 0 | 1 | 0 | 26 | 2 | 23 | 0 | 11 | 1 | 231 |
| toxic_channel | 76 | 57 | 0 | 0 | 0 | 59 | 89 | 0 | 0 | 7 | 0 | 0 | 0 | 1 | 0 | 23 | 1 | 18 | 0 | 7 | 1 | 220 |
| electron_yield | 37 | 40 | 0 | 0 | 0 | 22 | 42 | 0 | 8 | 1 | 11 | 10 | 8 | 3 | 6 | 22 | 1 | 2 | 0 | 6 | 8 | 179 |
| dispersal | 85 | 61 | 0 | 0 | 0 | 49 | 115 | 0 | 0 | 5 | 0 | 0 | 0 | 3 | 1 | 21 | 0 | 13 | 0 | 13 | 1 | 212 |
| forced | 67 | 50 | 0 | 0 | 0 | 51 | 78 | 0 | 0 | 5 | 0 | 0 | 0 | 3 | 1 | 18 | 1 | 17 | 0 | 8 | 1 | 149 |

## T8. Species-resolved chemical-potential budget (kJ per mmol C of NPP; residual relative to T*sigma_chem). T*sigma_chem = capture + (in - out) - reservoir - storage - other dissipation - residual

| config | n | T*sigma_chem / NPP | algal capture / NPP | boundary in - out / NPP | net to CO2/N2/H2O/H+ reservoirs / NPP | transport+hydrolysis+mortality+sediment / NPP | storage / NPP | residual / T*sigma_chem |
|---|---|---|---|---|---|---|---|---|
| emergent:neutral | 3 | 299.5 ± 5.9 | 508.8 ± 0.8 | 45.6 ± 7.8 | 248.2 ± 1.3 | 6.61 ± 0.30 | 0.10 ± 0.11 | 0.00011 ± 0.00010 |
| best_found_ref | 24 | 309.3 ± 6.7 | 511.1 ± 1.0 | 52.4 ± 6.1 | 250.6 ± 3.6 | 4.16 ± 0.56 | -0.03 ± 0.06 | -0.00180 ± 0.00119 |
| greedy | 24 | 296.3 ± 6.8 | 508.8 ± 0.7 | 43.9 ± 7.6 | 249.1 ± 1.4 | 7.12 ± 0.20 | 0.01 ± 0.01 | 0.00046 ± 0.00040 |
| random1 | 24 | 297.9 ± 10.2 | 511.0 ± 1.1 | 53.0 ± 4.9 | 261.0 ± 7.7 | 4.39 ± 0.38 | 0.03 ± 0.15 | 0.00215 ± 0.00301 |
| random2 | 24 | 299.7 ± 5.4 | 511.0 ± 0.8 | 52.5 ± 5.5 | 259.5 ± 5.4 | 4.35 ± 0.32 | -0.05 ± 0.12 | -0.00020 ± 0.00079 |
| emergent:crowd | 3 | 301.6 ± 4.2 | 508.5 ± 0.5 | 48.7 ± 5.7 | 249.1 ± 1.3 | 6.46 ± 0.32 | 0.03 ± 0.04 | 0.00007 ± 0.00007 |
| emergent:stackadv | 3 | 299.6 ± 6.1 | 508.8 ± 0.9 | 45.1 ± 8.2 | 247.5 ± 1.3 | 6.69 ± 0.33 | 0.10 ± 0.10 | 0.00012 ± 0.00011 |
| emergent:constitutive | 3 | 290.1 ± 2.0 | 510.2 ± 0.3 | 35.5 ± 2.7 | 249.2 ± 0.5 | 6.03 ± 0.14 | 0.28 ± 0.03 | 0.00028 ± 0.00006 |
| emergent:channel | 3 | 299.9 ± 5.7 | 508.8 ± 0.8 | 45.7 ± 7.7 | 248.0 ± 1.4 | 6.50 ± 0.38 | 0.08 ± 0.09 | 0.00011 ± 0.00010 |
| emergent:toxic_channel | 3 | 300.5 ± 5.0 | 508.7 ± 0.7 | 46.5 ± 7.0 | 248.1 ± 1.4 | 6.48 ± 0.41 | 0.06 ± 0.07 | 0.00010 ± 0.00009 |
| emergent:electron_yield | 3 | 303.9 ± 3.0 | 508.8 ± 0.6 | 51.9 ± 4.8 | 250.4 ± 1.9 | 6.36 ± 0.32 | -0.01 ± 0.09 | 0.00003 ± 0.00005 |
| emergent:dispersal | 3 | 300.7 ± 5.1 | 508.8 ± 0.5 | 48.7 ± 5.7 | 250.1 ± 0.1 | 6.55 ± 0.42 | 0.05 ± 0.07 | 0.00008 ± 0.00008 |
| emergent:forced | 3 | 290.4 ± 12.7 | 507.7 ± 0.6 | 35.1 ± 10.8 | 262.5 ± 8.4 | 5.85 ± 0.58 | -16.03 ± 5.39 | 0.00010 ± 0.00009 |
