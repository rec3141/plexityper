"""Condense results/*.tsv into results/summary_tables.md (markdown tables used in README)."""
import os, sys, csv
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, sys.argv[1] if len(sys.argv) > 1 else 'results')


def read(name):
    p = os.path.join(RES, name + '.tsv')
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


LABELS = {  # revised display labels (README section 6; Astra 249011970736832512)
    'mean_n_required': 'mean distinct donors in witnessed supports (not partners needed)',
    'kappa_all_supports_frac': 'certified all-support cycle frac (none certifiable: supports incomplete)',
    'kappa_any_support_frac': 'witnessed cyclic support frac',
}


def md(rows, keys, fmt=None):
    if not rows:
        return '_(no rows)_\n'
    out = '| ' + ' | '.join(LABELS.get(k, k) for k in keys) + ' |\n|' + '---|' * len(keys) + '\n'
    for r in rows:
        vals = []
        for k in keys:
            v = r.get(k, '')
            try:
                f = float(v)
                v = ('%.3g' % f) if abs(f) < 1e4 else ('%.2e' % f)
            except (ValueError, TypeError):
                pass
            vals.append(str(v))
        out += '| ' + ' | '.join(vals) + ' |\n'
    return out


def main():
    S = []
    S.append('## Truth summary (viable lineage-times among the 300 most abundant, every 14 d)\n')
    tr = read('truth_summary')
    keys = ['world', 'n_viable', 'genome_dep_d_h1>=1_frac', 'plexi_regulated_d_h1>=2_frac',
            'd>=2_functional_fan_in_frac', 'd>=2_dilution_witnessed_frac', 'd>=2_dilution_unresolved_frac',
            'd>=2_other_frac', 'd_hi_unresolved_frac', 'kappa_any_support_frac', 'kappa_all_supports_frac',
            'layer_expression_frac', 'euk_required_frac', 'realised_not_required>=1_frac', 'mean_n_realised',
            'mean_n_required', 'rows_with_d>=1', 'supports_unknown_rows']
    S.append(md(tr, keys))
    S.append('\nNotes: no certified all-support cycle (all d >= 2 searches are incomplete, so the '
             'certified field cannot be 1); witnessed cyclic supports in 4 W_powerlaw / 11 W_lognormal rows. '
             '"mean distinct donors in witnessed supports" averages the union of witnessed donors over all '
             'viable rows (zeros for unknown supports; interchangeable donors counted separately) and does '
             'not estimate how many partners a recipient needs.\n')
    S.append('\n### Genes lost and Liebig masking\n')
    S.append(md(tr, ['world', 'genes_lost_0_frac', 'd>=2_given_genes_lost_0', 'genes_lost_1_frac',
                     'd>=2_given_genes_lost_1', 'genes_lost_2-4_frac', 'd>=2_given_genes_lost_2-4',
                     'genes_lost_5_frac', 'd>=2_given_genes_lost_5', 'd>=2_with_m0_deficient_frac']))
    S.append('\n## Knock-out validation (dynamic 20-d removal vs static call)\n')
    S.append(md(read('knockout_validation_summary'), ['world', 'test', 'static_class', 'n_pairs', 'median_ko_ratio',
                                                      'frac_ratio_below_0_5', 'frac_ratio_below_0_8']))
    S.append('\n## Timescales by design and layer\n')
    ts = read('timescales_by_design_layer')
    S.append(md(ts, ['world', 'design', 'layer', 'n_series', 'n_times', 'median_dt_d', 'tau_hat_median',
                     'tau_hat_q25', 'tau_hat_q75', 'frac_tau_at_bounds', 'ar1_null_reject_frac',
                     'power_fast_lt30d', 'power_mid_30_180d', 'power_slow_gt180d', 'auc_vs_next_faster']))
    S.append('\n## Eukaryote bloom dominance (daily truth)\n')
    S.append(md(read('eukaryote_bloom_dominance'), ['world', 'layer', 'n_series', 'r2_season_T',
                                                   'r2_gain_lagged_algae', 'r2_gain_q75', 'power_fast',
                                                   'power_mid', 'power_slow', 'algal_share_of_organic_C_supply',
                                                   'algal_biomass_cv']))
    S.append('\n## Perturbation propagation (two campaigns)\n')
    S.append(md(read('perturbation_propagation'), ['campaign', 'source', 'layer', 't10_d', 't50_d', 'tpeak_d',
                                                  'peak_div', 'pre_div_mean', 'pre_div_sd']))
    S.append('\n## Aliasing\n')
    S.append(md(read('aliasing'), ['world', 'design', 'normalisation', 'n', 'peak_period_d', 'power_frac_28_42d',
                                   'r2_by_sampling_hour', 'tau_hat_median']))
    S.append('\n## Closure\n')
    S.append(md(read('closure'), ['world', 'measure', 'n_taxa', 'median_r_rel_vs_abs', 'mean_pairwise_corr_abs',
                                  'mean_pairwise_corr_rel', 'frac_pairs_sign_flip', 'frac_pairs_spurious_neg',
                                  'median_r_with_euk_read_frac_rel', 'median_r_with_euk_read_frac_abs',
                                  'slow_power_abs', 'slow_power_rel']))
    S.append('\n## P3 pilot (protocol-aligned where feasible): verdict vs oracle verdict\n')
    p3 = read('p3_groups')
    agg = defaultdict(Counter)
    for r in p3:
        agg[(r['world'], r['design'], r['universe'])][(r.get('verdict', ''), r.get('oracle_verdict', ''))] += 1
    rows = []
    for k in sorted(agg):
        for (v, o), n in sorted(agg[k].items()):
            rows.append(dict(world=k[0], design=k[1], universe=k[2], pilot_verdict=v, oracle_verdict=o, n=n))
    S.append(md(rows, ['world', 'design', 'universe', 'pilot_verdict', 'oracle_verdict', 'n']))
    S.append('\n### P3 groups (all)\n')
    S.append(md(p3, ['world', 'design', 'universe', 'n_members', 'verdict', 'gain_me', 'gain_me_lo', 'gain_me_hi',
                     'gain_mn', 'gain_mn_lo', 'gain_mn_hi', 'oracle_verdict', 'oracle_gain_me', 'oracle_gain_me_lo',
                     'oracle_gain_mn', 'oracle_gain_mn_lo', 'truth_direct_required_edge_within_group', 'n_test']))
    S.append('\n## Abundant core (frozen on training events, equal event weights)\n')
    S.append(md(read('core_ranking'), ['world', 'design', 'n_train', 'core_n', 'core_read_share',
                                       'core_detect_freq_median', 'core_seasonal_occ_median',
                                       'core_frac_occ_lt_0_5', 'noncore_detected_taxa']))
    S.append('\n## Bloomers vs core: genome-layer dependencies in truth\n')
    S.append(md(read('bloomers_vs_core'), ['world', 'group', 'n', 'n_with_truth', 'frac_with_genome_dependency',
                                           'frac_dependency_and_genes_lost_ge1', 'frac_functional_fan_in']))
    S.append('\n## P1 song vs singer\n')
    S.append(md(read('p1_song_vs_singer'), ['study', 'world', 'guild_set', 'n_stations', 'n_obs', 'gate_pass',
                                            'n_dominant_genera_ge5_stations', 'transfer_genera_tested',
                                            'transfer_genera_gain_lo_gt_margin', 'guild_gain',
                                            'guild_gain_lo', 'guild_gain_hi', 'identity_gain', 'identity_gain_lo',
                                            'identity_gain_hi', 'guild_slope_bootstrap_hi',
                                            'role_recurrence_decision', 'identity_decision',
                                            'truth_identity_matters', 'truth_sd_usable_fraction',
                                            'truth_r_pool_vs_uptake_expr', 'r_guild_vs_true_pool',
                                            'r_observed_y_vs_true_uptake', 'sd_ratio_observed_to_true']))
    S.append('\n## Mechinterp probes\n')
    S.append(md(read('mechinterp_probes'), ['world', 'probe', 'target', 'kind', 'score', 'null_mean',
                                            'percentile_vs_random', 'size']))
    open(os.path.join(RES, 'summary_tables.md'), 'w').write('\n'.join(S))
    print('wrote', os.path.join(RES, 'summary_tables.md'))


if __name__ == '__main__':
    main()
