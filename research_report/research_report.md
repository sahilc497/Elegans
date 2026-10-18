# C. elegans Connectome Simulation: Results Report

## 1. Scope

This report summarizes the existing simulation output files. No simulations are rerun by this script.

The model uses the White et al. 1986 connectome loaded through CECT. It is an exploratory rate-based network model, not a validated physiological model of the worm.

## 2. Robustness analysis

- Parameter combinations represented: 27
- Mean left/right cosine similarity: 0.388997
- Minimum similarity: 0.032733
- Maximum similarity: 0.995744

### By chemical strength

|   chemical_strength |   count |      mean |       min |       max |
|--------------------:|--------:|----------:|----------:|----------:|
|                 0.4 |       9 | 0.0627238 | 0.0327332 | 0.0956865 |
|                 0.8 |       9 | 0.167377  | 0.0774338 | 0.299902  |
|                 1.2 |       9 | 0.936889  | 0.788042  | 0.995744  |

## 3. Coupling ablation

The following table compares the full model with chemical-only, gap-only, and no-recurrent-coupling configurations.

| mode                  |   combinations |   mean_similarity |   min_similarity |   max_similarity |   mean_left_integrated |   mean_right_integrated |
|:----------------------|---------------:|------------------:|-----------------:|-----------------:|-----------------------:|------------------------:|
| full                  |             27 |         0.389036  |       0.0327332  |        0.995761  |              15.7709   |               17.0499   |
| chemical_only         |             27 |         0.290778  |       0.00254965 |        0.960723  |               7.91665  |                7.92285  |
| gap_only              |             27 |         0.0552668 |       0.0287901  |        0.0831036 |               0.987505 |                1.04707  |
| no_recurrent_coupling |             27 |         0         |       0          |        0         |               0.904043 |                0.904043 |

## 4. Propagation summary

| stimulus   | mode          |   peak_network_activity |   peak_mean_activity |   integrated_network_activity |   maximum_active_neurons |   maximum_excess_active_neurons |   neurons_above_threshold_at_end |
|:-----------|:--------------|------------------------:|---------------------:|------------------------------:|-------------------------:|--------------------------------:|---------------------------------:|
| ASHL       | full          |                0.676351 |           0.00521233 |                      3.93476  |                      258 |                             257 |                                6 |
| ASHL       | chemical_only |                0.784328 |           0.00475039 |                      2.75337  |                      200 |                             200 |                                2 |
| ASHL       | gap_only      |                0.6748   |           0.00266125 |                      0.987235 |                       26 |                              25 |                                0 |
| ASHR       | full          |                0.677462 |           0.00558149 |                      4.31199  |                      259 |                             258 |                                8 |
| ASHR       | chemical_only |                0.783814 |           0.00474559 |                      2.74082  |                      195 |                             194 |                                4 |
| ASHR       | gap_only      |                0.675594 |           0.00280513 |                      1.04879  |                       25 |                              24 |                                0 |

## 5. Threshold sensitivity

The number of neurons counted as active changes with the chosen numerical threshold. These thresholds are analytical cutoffs, not experimentally established firing thresholds.

| stimulus   | mode          |   threshold |   neurons_ever_above_threshold |   neurons_above_threshold_at_end |   neurons_with_integrated_activity_above_threshold |
|:-----------|:--------------|------------:|-------------------------------:|---------------------------------:|---------------------------------------------------:|
| ASHL       | full          |      1e-06  |                            295 |                              229 |                                                295 |
| ASHL       | full          |      1e-05  |                            289 |                               87 |                                                293 |
| ASHL       | full          |      0.0001 |                            258 |                                6 |                                                272 |
| ASHL       | full          |      0.001  |                            181 |                                0 |                                                234 |
| ASHL       | full          |      0.01   |                             25 |                                0 |                                                 96 |
| ASHL       | full          |      0.05   |                              4 |                                0 |                                                 11 |
| ASHL       | chemical_only |      1e-06  |                            264 |                               69 |                                                266 |
| ASHL       | chemical_only |      1e-05  |                            248 |                               19 |                                                259 |
| ASHL       | chemical_only |      0.0001 |                            201 |                                2 |                                                228 |
| ASHL       | chemical_only |      0.001  |                            108 |                                0 |                                                148 |
| ASHL       | chemical_only |      0.01   |                             16 |                                0 |                                                 39 |
| ASHL       | chemical_only |      0.05   |                              5 |                                0 |                                                  9 |
| ASHL       | gap_only      |      1e-06  |                             70 |                                0 |                                                 81 |
| ASHL       | gap_only      |      1e-05  |                             50 |                                0 |                                                 54 |
| ASHL       | gap_only      |      0.0001 |                             26 |                                0 |                                                 31 |
| ASHL       | gap_only      |      0.001  |                             14 |                                0 |                                                 17 |
| ASHL       | gap_only      |      0.01   |                              6 |                                0 |                                                  6 |
| ASHL       | gap_only      |      0.05   |                              1 |                                0 |                                                  2 |
| ASHR       | full          |      1e-06  |                            293 |                              235 |                                                295 |
| ASHR       | full          |      1e-05  |                            276 |                               93 |                                                288 |
| ASHR       | full          |      0.0001 |                            259 |                                8 |                                                270 |
| ASHR       | full          |      0.001  |                            185 |                                0 |                                                236 |
| ASHR       | full          |      0.01   |                             34 |                                0 |                                                103 |
| ASHR       | full          |      0.05   |                              4 |                                0 |                                                 12 |
| ASHR       | chemical_only |      1e-06  |                            260 |                               67 |                                                266 |
| ASHR       | chemical_only |      1e-05  |                            236 |                               16 |                                                253 |
| ASHR       | chemical_only |      0.0001 |                            195 |                                4 |                                                219 |
| ASHR       | chemical_only |      0.001  |                            107 |                                0 |                                                144 |
| ASHR       | chemical_only |      0.01   |                             16 |                                0 |                                                 47 |
| ASHR       | chemical_only |      0.05   |                              5 |                                0 |                                                  8 |
| ASHR       | gap_only      |      1e-06  |                             70 |                                0 |                                                 82 |
| ASHR       | gap_only      |      1e-05  |                             50 |                                0 |                                                 57 |
| ASHR       | gap_only      |      0.0001 |                             25 |                                0 |                                                 31 |
| ASHR       | gap_only      |      0.001  |                             15 |                                0 |                                                 19 |
| ASHR       | gap_only      |      0.01   |                              7 |                                0 |                                                  7 |
| ASHR       | gap_only      |      0.05   |                              2 |                                0 |                                                  2 |

## 6. Figures

- `figure_robustness.png`: response similarity by chemical strength.
- `figure_coupling_ablation.png`: response similarity across coupling modes.
- `figure_threshold_sensitivity.png`: active-neuron count versus threshold.

## 7. Interpretation and limitations

1. Chemical-only coupling produced broader activity propagation than gap-only coupling under the tested equations and parameters.
2. Left/right stimulation produced different neuron-level response profiles in this implementation.
3. Similarity and active-neuron counts depend on model parameters and the activity threshold.
4. The model treats chemical connections as positive; real chemical synapses cannot all be assumed excitatory.
5. Gap-junction symmetrization, connection normalization, activity dynamics, and parameter values are modeling assumptions.
6. These results establish behavior of this implementation, not experimental biological causality or validated worm behavior.

## 8. Reproducibility

The CSV tables in this folder are derived from the existing result files. The original simulation scripts and source CSVs should be retained alongside this report for full reproducibility.
