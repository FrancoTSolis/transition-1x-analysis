
### smallval: protocol lin, reference fci

| candidate | n | LUCJ % CCSD corr | LUCJ err (mHa) | QSCI err seed 0 | QSCI err 5-seed | seed sd | QSCI % corr | dim/spin | p_HF | configs in 10^5 | QSCI better than label |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 9 | 16.6 | 197.3 | 82.37 | 85.02 | 7.96 | 63.62 | 83 | 0.968 | 173 | 0/9 |
| compressed DF (Lin) | 9 | -401.3 | 1127.5 | 8.00 | 8.00 | 0.25 | 96.61 | 929 | 0.320 | 10433 | 7/9 |
| optimize=True label | 9 | 62.0 | 94.7 | 9.81 | 10.10 | 0.46 | 95.70 | 520 | 0.842 | 2097 | — |
| pretrained x4 | 9 | 59.4 | 99.5 | 11.44 | 11.41 | 0.41 | 95.14 | 448 | 0.854 | 1809 | 1/9 |
| NN+RL rl4L | 9 | 68.2 | 80.1 | 8.41 | 8.61 | 0.26 | 96.32 | 611 | 0.839 | 2391 | 8/9 |
| NN+RL rl4n29f | 9 | 70.3 | 75.5 | 9.93 | 9.97 | 0.45 | 95.74 | 538 | 0.860 | 2137 | 4/9 |

QSCI error ratio truncated / candidate (5-seed means, reference fci):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| compressed DF (Lin) | 9 | 12.4 | 11.2 | 4.8 | 22.6 |
| optimize=True label | 9 | 8.6 | 8.4 | 5.1 | 12.2 |
| pretrained x4 | 9 | 7.7 | 7.4 | 4.3 | 10.7 |
| NN+RL rl4L | 9 | 10.5 | 9.9 | 5.5 | 15.2 |
| NN+RL rl4n29f | 9 | 9.1 | 8.6 | 4.6 | 14.0 |

QSCI error ratio optimize=True label / candidate (reference fci; > 1 = candidate better; molecules with an error <= 0 skipped):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| compressed DF (Lin) | 9 | 1.39 | 1.34 | 0.81 | 1.90 |
| pretrained x4 | 9 | 0.89 | 0.89 | 0.81 | 1.00 |
| NN+RL rl4L | 9 | 1.21 | 1.19 | 0.86 | 1.61 |
| NN+RL rl4n29f | 9 | 1.04 | 1.03 | 0.78 | 1.29 |

QSCI error ratio compressed DF (Lin) / candidate (reference fci; > 1 = candidate better; molecules with an error <= 0 skipped):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| optimize=True label | 9 | 0.78 | 0.75 | 0.53 | 1.23 |
| pretrained x4 | 9 | 0.70 | 0.67 | 0.44 | 1.23 |
| NN+RL rl4L | 9 | 0.94 | 0.89 | 0.57 | 1.64 |
| NN+RL rl4n29f | 9 | 0.82 | 0.77 | 0.51 | 1.43 |

### smallval: protocol lin, reference ccsdt

| candidate | n | LUCJ % CCSD corr | LUCJ err (mHa) | QSCI err seed 0 | QSCI err 5-seed | seed sd | QSCI % corr | dim/spin | p_HF | configs in 10^5 | QSCI better than label |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 9 | 16.6 | 195.4 | 80.46 | 83.11 | 7.96 | 64.13 | 83 | 0.968 | 173 | 0/9 |
| compressed DF (Lin) | 9 | -401.3 | 1125.6 | 6.09 | 6.09 | 0.25 | 97.38 | 929 | 0.320 | 10433 | 7/9 |
| optimize=True label | 9 | 62.0 | 92.8 | 7.90 | 8.19 | 0.46 | 96.46 | 520 | 0.842 | 2097 | — |
| pretrained x4 | 9 | 59.4 | 97.6 | 9.53 | 9.50 | 0.41 | 95.90 | 448 | 0.854 | 1809 | 1/9 |
| NN+RL rl4L | 9 | 68.2 | 78.2 | 6.50 | 6.70 | 0.26 | 97.09 | 611 | 0.839 | 2391 | 8/9 |
| NN+RL rl4n29f | 9 | 70.3 | 73.6 | 8.02 | 8.06 | 0.45 | 96.50 | 538 | 0.860 | 2137 | 4/9 |

QSCI error ratio truncated / candidate (5-seed means, reference ccsdt):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| compressed DF (Lin) | 9 | 18.6 | 15.3 | 5.4 | 43.8 |
| optimize=True label | 9 | 10.5 | 10.1 | 5.8 | 14.7 |
| pretrained x4 | 9 | 9.0 | 8.7 | 4.7 | 11.6 |
| NN+RL rl4L | 9 | 13.7 | 12.8 | 6.3 | 19.4 |
| NN+RL rl4n29f | 9 | 11.3 | 10.5 | 5.1 | 15.7 |

QSCI error ratio optimize=True label / candidate (reference ccsdt; > 1 = candidate better; molecules with an error <= 0 skipped):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| compressed DF (Lin) | 9 | 1.65 | 1.51 | 0.79 | 2.98 |
| pretrained x4 | 9 | 0.87 | 0.86 | 0.77 | 1.00 |
| NN+RL rl4L | 9 | 1.29 | 1.26 | 0.84 | 1.85 |
| NN+RL rl4n29f | 9 | 1.06 | 1.04 | 0.75 | 1.38 |

QSCI error ratio compressed DF (Lin) / candidate (reference ccsdt; > 1 = candidate better; molecules with an error <= 0 skipped):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| optimize=True label | 9 | 0.72 | 0.66 | 0.34 | 1.27 |
| pretrained x4 | 9 | 0.63 | 0.57 | 0.26 | 1.28 |
| NN+RL rl4L | 9 | 0.93 | 0.83 | 0.44 | 1.81 |
| NN+RL rl4n29f | 9 | 0.77 | 0.69 | 0.33 | 1.52 |

### smallval: protocol lin, candidate vs label (QSCI: 5-seed means)

| candidate vs optimize=True label | n | mean dE_LUCJ (mHa) | LUCJ lower | mean dE_QSCI (mHa) | QSCI lower | QSCI lower by > 2 s.e. | QSCI higher by > 2 s.e. |
|:--|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 9 | +102.6 | 0/9 | +74.92 | 0/9 | 0/9 | 9/9 |
| compressed DF (Lin) | 9 | +1032.7 | 0/9 | -2.09 | 7/9 | 7/9 | 2/9 |
| pretrained x4 | 9 | +4.8 | 5/9 | +1.32 | 1/9 | 0/9 | 6/9 |
| NN+RL rl4L | 9 | -14.6 | 8/9 | -1.48 | 8/9 | 8/9 | 1/9 |
| NN+RL rl4n29f | 9 | -19.2 | 9/9 | -0.12 | 4/9 | 4/9 | 3/9 |

* ranking LUCJ vs QSCI (all 6): rho per molecule mean 0.13 (median 0.03, min -0.26, max 0.89), pooled -0.09; concordant pairs 77/135 (resolved beyond 2 s.e.: 75/128); per molecule [0.89, -0.14, -0.26, -0.14, -0.03, 0.2, 0.03, 0.09, 0.54]
  pooled Spearman of QSCI energy with dim/spin -0.638, p_HF 0.665, entropy -0.646, LUCJ energy -0.09 (n=54)
* ranking LUCJ vs QSCI (without cdf_lin): rho per molecule mean 0.66 (median 0.70, min 0.30, max 0.90), pooled 0.48; concordant pairs 69/90 (resolved beyond 2 s.e.: 67/85); per molecule [0.9, 0.5, 0.3, 0.5, 0.7, 0.5, 0.8, 0.9, 0.8]
  pooled Spearman of QSCI energy with dim/spin -0.68, p_HF 0.64, entropy -0.651, LUCJ energy 0.478 (n=45)
* ranking LUCJ vs QSCI (label + network family): rho per molecule mean 0.31 (median 0.40, min -0.40, max 0.80), pooled 0.24; concordant pairs 33/54 (resolved beyond 2 s.e.: 31/49); per molecule [0.8, 0.0, -0.4, 0.0, 0.4, 0.0, 0.6, 0.8, 0.6]
  pooled Spearman of QSCI energy with dim/spin -0.957, p_HF 0.685, entropy -0.872, LUCJ energy 0.241 (n=36)

### smallval: protocol n2631g, reference fci

| candidate | n | LUCJ % CCSD corr | LUCJ err (mHa) | QSCI err seed 0 | QSCI err 5-seed | seed sd | QSCI % corr | dim/spin | p_HF | configs in 10^5 | QSCI better than label |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 9 | 16.6 | 197.3 | 48.38 | 48.15 | 2.22 | 79.38 | 174 | 0.968 | 173 | 0/9 |
| compressed DF (Lin) | 9 | -401.3 | 1127.5 | 0.69 | 0.69 | – | 99.71 | 3277 | 0.320 | 10433 | 9/9 |
| optimize=True label | 9 | 62.0 | 94.7 | 2.51 | 2.55 | 0.10 | 98.91 | 1192 | 0.842 | 2097 | — |
| pretrained x4 | 9 | 59.4 | 99.5 | 3.12 | 3.12 | 0.12 | 98.67 | 1067 | 0.854 | 1809 | 0/9 |
| NN+RL rl4L | 9 | 68.2 | 80.1 | 2.05 | 2.08 | 0.06 | 99.11 | 1407 | 0.839 | 2391 | 8/9 |
| NN+RL rl4n29f | 9 | 70.3 | 75.5 | 2.52 | 2.58 | 0.08 | 98.90 | 1242 | 0.860 | 2137 | 6/9 |

QSCI error ratio truncated / candidate (5-seed means, reference fci):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| compressed DF (Lin) | 9 | 89.6 | 77.3 | 40.7 | 228.7 |
| optimize=True label | 9 | 20.1 | 18.1 | 7.5 | 40.0 |
| pretrained x4 | 9 | 16.3 | 14.8 | 6.1 | 29.3 |
| NN+RL rl4L | 9 | 26.0 | 22.8 | 8.8 | 53.4 |
| NN+RL rl4n29f | 9 | 21.0 | 18.3 | 6.9 | 49.4 |

QSCI error ratio optimize=True label / candidate (reference fci; > 1 = candidate better; molecules with an error <= 0 skipped):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| compressed DF (Lin) | 9 | 5.00 | 4.27 | 1.64 | 10.57 |
| pretrained x4 | 9 | 0.82 | 0.82 | 0.73 | 0.91 |
| NN+RL rl4L | 9 | 1.28 | 1.26 | 0.84 | 1.68 |
| NN+RL rl4n29f | 9 | 1.03 | 1.01 | 0.73 | 1.25 |

QSCI error ratio compressed DF (Lin) / candidate (reference fci; > 1 = candidate better; molecules with an error <= 0 skipped):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| optimize=True label | 9 | 0.28 | 0.23 | 0.09 | 0.61 |
| pretrained x4 | 9 | 0.23 | 0.19 | 0.07 | 0.45 |
| NN+RL rl4L | 9 | 0.36 | 0.30 | 0.12 | 0.81 |
| NN+RL rl4n29f | 9 | 0.30 | 0.24 | 0.10 | 0.75 |

### smallval: protocol n2631g, reference ccsdt

| candidate | n | LUCJ % CCSD corr | LUCJ err (mHa) | QSCI err seed 0 | QSCI err 5-seed | seed sd | QSCI % corr | dim/spin | p_HF | configs in 10^5 | QSCI better than label |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 9 | 16.6 | 195.4 | 46.47 | 46.24 | 2.22 | 80.02 | 174 | 0.968 | 173 | 0/9 |
| compressed DF (Lin) | 9 | -401.3 | 1125.6 | -1.22 | -1.22 | – | 100.50 | 3277 | 0.320 | 10433 | 9/9 |
| optimize=True label | 9 | 62.0 | 92.8 | 0.60 | 0.64 | 0.10 | 99.70 | 1192 | 0.842 | 2097 | — |
| pretrained x4 | 9 | 59.4 | 97.6 | 1.21 | 1.21 | 0.12 | 99.46 | 1067 | 0.854 | 1809 | 0/9 |
| NN+RL rl4L | 9 | 68.2 | 78.2 | 0.14 | 0.17 | 0.06 | 99.90 | 1407 | 0.839 | 2391 | 8/9 |
| NN+RL rl4n29f | 9 | 70.3 | 73.6 | 0.61 | 0.67 | 0.08 | 99.68 | 1242 | 0.860 | 2137 | 6/9 |

QSCI error ratio truncated / candidate (5-seed means, reference ccsdt):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| compressed DF (Lin) | 1 | 197.6 | 197.6 | 197.6 | 197.6 |
| optimize=True label | 7 | 48.4 | 43.7 | 16.7 | 69.3 |
| pretrained x4 | 8 | 37.3 | 32.3 | 10.7 | 73.4 |
| NN+RL rl4L | 6 | 153.4 | 74.0 | 21.1 | 628.5 |
| NN+RL rl4n29f | 7 | 76.6 | 47.3 | 13.6 | 261.0 |

QSCI error ratio optimize=True label / candidate (reference ccsdt; > 1 = candidate better; molecules with an error <= 0 skipped):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| compressed DF (Lin) | 1 | 2.94 | 2.94 | 2.94 | 2.94 |
| pretrained x4 | 7 | 0.66 | 0.66 | 0.61 | 0.74 |
| NN+RL rl4L | 6 | 2.81 | 1.83 | 0.69 | 10.01 |
| NN+RL rl4n29f | 7 | 1.35 | 1.08 | 0.55 | 3.76 |

QSCI error ratio compressed DF (Lin) / candidate (reference ccsdt; > 1 = candidate better; molecules with an error <= 0 skipped):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| optimize=True label | 1 | 0.34 | 0.34 | 0.34 | 0.34 |
| pretrained x4 | 1 | 0.21 | 0.21 | 0.21 | 0.21 |
| NN+RL rl4L | 1 | 0.59 | 0.59 | 0.59 | 0.59 |
| NN+RL rl4n29f | 1 | 0.50 | 0.50 | 0.50 | 0.50 |

### smallval: protocol n2631g, candidate vs label (QSCI: 5-seed means)

| candidate vs optimize=True label | n | mean dE_LUCJ (mHa) | LUCJ lower | mean dE_QSCI (mHa) | QSCI lower | QSCI lower by > 2 s.e. | QSCI higher by > 2 s.e. |
|:--|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 9 | +102.6 | 0/9 | +45.60 | 0/9 | 0/9 | 9/9 |
| compressed DF (Lin) | 9 | +1032.7 | 0/9 | -1.86 | 9/9 | 9/9 | 0/9 |
| pretrained x4 | 9 | +4.8 | 5/9 | +0.57 | 0/9 | 0/9 | 9/9 |
| NN+RL rl4L | 9 | -14.6 | 8/9 | -0.47 | 8/9 | 8/9 | 1/9 |
| NN+RL rl4n29f | 9 | -19.2 | 9/9 | +0.03 | 6/9 | 3/9 | 3/9 |

* ranking LUCJ vs QSCI (all 6): rho per molecule mean -0.03 (median 0.03, min -0.26, max 0.09), pooled -0.30; concordant pairs 70/135 (resolved beyond 2 s.e.: 67/132); per molecule [0.03, 0.03, -0.26, -0.14, 0.09, -0.14, 0.03, 0.09, 0.03]
  pooled Spearman of QSCI energy with dim/spin -0.808, p_HF 0.788, entropy -0.812, LUCJ energy -0.298 (n=54)
* ranking LUCJ vs QSCI (without cdf_lin): rho per molecule mean 0.70 (median 0.80, min 0.30, max 0.90), pooled 0.49; concordant pairs 70/90 (resolved beyond 2 s.e.: 67/87); per molecule [0.8, 0.8, 0.3, 0.5, 0.9, 0.5, 0.8, 0.9, 0.8]
  pooled Spearman of QSCI energy with dim/spin -0.579, p_HF 0.548, entropy -0.546, LUCJ energy 0.492 (n=45)
* ranking LUCJ vs QSCI (label + network family): rho per molecule mean 0.40 (median 0.60, min -0.40, max 0.80), pooled 0.29; concordant pairs 34/54 (resolved beyond 2 s.e.: 31/51); per molecule [0.6, 0.6, -0.4, 0.0, 0.8, 0.0, 0.6, 0.8, 0.6]
  pooled Spearman of QSCI energy with dim/spin -0.949, p_HF 0.634, entropy -0.842, LUCJ energy 0.286 (n=36)

### smallval: protocol top100, reference fci

| candidate | n | LUCJ % CCSD corr | LUCJ err (mHa) | QSCI err seed 0 | QSCI err 5-seed | seed sd | QSCI % corr | dim/spin | p_HF | configs in 10^5 | QSCI better than label |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 9 | 16.6 | 197.3 | 64.49 | 64.49 | – | 72.30 | 100 | 0.968 | 173 | 0/9 |
| compressed DF (Lin) | 9 | -401.3 | 1127.5 | 37.63 | 37.63 | – | 83.99 | 100 | 0.320 | 10433 | 1/9 |
| optimize=True label | 9 | 62.0 | 94.7 | 29.57 | 29.57 | – | 87.49 | 100 | 0.842 | 2097 | — |
| pretrained x4 | 9 | 59.4 | 99.5 | 28.96 | 28.96 | – | 87.74 | 100 | 0.854 | 1809 | 6/9 |
| NN+RL rl4L | 9 | 68.2 | 80.1 | 29.80 | 29.80 | – | 87.36 | 100 | 0.839 | 2391 | 5/9 |
| NN+RL rl4n29f | 9 | 70.3 | 75.5 | 30.57 | 30.57 | – | 87.03 | 100 | 0.860 | 2137 | 2/9 |

### smallval: protocol top100, reference ccsdt

| candidate | n | LUCJ % CCSD corr | LUCJ err (mHa) | QSCI err seed 0 | QSCI err 5-seed | seed sd | QSCI % corr | dim/spin | p_HF | configs in 10^5 | QSCI better than label |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 9 | 16.6 | 195.4 | 62.58 | 62.58 | – | 72.89 | 100 | 0.968 | 173 | 0/9 |
| compressed DF (Lin) | 9 | -401.3 | 1125.6 | 35.72 | 35.72 | – | 84.66 | 100 | 0.320 | 10433 | 1/9 |
| optimize=True label | 9 | 62.0 | 92.8 | 27.66 | 27.66 | – | 88.18 | 100 | 0.842 | 2097 | — |
| pretrained x4 | 9 | 59.4 | 97.6 | 27.05 | 27.05 | – | 88.44 | 100 | 0.854 | 1809 | 6/9 |
| NN+RL rl4L | 9 | 68.2 | 78.2 | 27.89 | 27.89 | – | 88.06 | 100 | 0.839 | 2391 | 5/9 |
| NN+RL rl4n29f | 9 | 70.3 | 73.6 | 28.66 | 28.66 | – | 87.72 | 100 | 0.860 | 2137 | 2/9 |

### smallval: protocol top100, candidate vs label (QSCI: 5-seed means)

| candidate vs optimize=True label | n | mean dE_LUCJ (mHa) | LUCJ lower | mean dE_QSCI (mHa) | QSCI lower | QSCI lower by > 2 s.e. | QSCI higher by > 2 s.e. |
|:--|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 9 | +102.6 | 0/9 | +34.92 | 0/9 | 0/9 | 9/9 |
| compressed DF (Lin) | 9 | +1032.7 | 0/9 | +8.05 | 1/9 | 1/9 | 8/9 |
| pretrained x4 | 9 | +4.8 | 5/9 | -0.61 | 6/9 | 6/9 | 3/9 |
| NN+RL rl4L | 9 | -14.6 | 8/9 | +0.22 | 5/9 | 5/9 | 4/9 |
| NN+RL rl4n29f | 9 | -19.2 | 9/9 | +1.00 | 2/9 | 2/9 | 7/9 |

* ranking LUCJ vs QSCI (all 6): rho per molecule mean 0.47 (median 0.49, min 0.09, max 0.77), pooled 0.51; concordant pairs 89/135 (resolved beyond 2 s.e.: 89/135); per molecule [0.66, 0.2, 0.09, 0.49, 0.37, 0.54, 0.77, 0.77, 0.37]
  pooled Spearman of QSCI energy with dim/spin -0.43, p_HF 0.099, entropy -0.129, LUCJ energy 0.511 (n=54)
* ranking LUCJ vs QSCI (without cdf_lin): rho per molecule mean 0.33 (median 0.30, min 0.00, max 0.70), pooled 0.41; concordant pairs 56/90 (resolved beyond 2 s.e.: 56/90); per molecule [0.4, 0.4, 0.3, 0.2, 0.0, 0.3, 0.7, 0.7, 0.0]
  pooled Spearman of QSCI energy with dim/spin -0.545, p_HF 0.489, entropy -0.472, LUCJ energy 0.412 (n=45)
* ranking LUCJ vs QSCI (label + network family): rho per molecule mean -0.33 (median -0.40, min -1.00, max 0.40), pooled -0.33; concordant pairs 20/54 (resolved beyond 2 s.e.: 20/54); per molecule [-0.2, -0.2, -0.4, -0.6, -1.0, -0.4, 0.4, 0.4, -1.0]
  pooled Spearman of QSCI energy with dim/spin nan, p_HF 0.162, entropy -0.068, LUCJ energy -0.335 (n=36)

### smallval: protocol top300, reference fci

| candidate | n | LUCJ % CCSD corr | LUCJ err (mHa) | QSCI err seed 0 | QSCI err 5-seed | seed sd | QSCI % corr | dim/spin | p_HF | configs in 10^5 | QSCI better than label |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 9 | 16.6 | 197.3 | 49.38 | 49.38 | – | 78.89 | 167 | 0.968 | 173 | 0/9 |
| compressed DF (Lin) | 9 | -401.3 | 1127.5 | 18.12 | 18.12 | – | 92.33 | 300 | 0.320 | 10433 | 0/9 |
| optimize=True label | 9 | 62.0 | 94.7 | 13.89 | 13.89 | – | 94.11 | 300 | 0.842 | 2097 | — |
| pretrained x4 | 9 | 59.4 | 99.5 | 13.43 | 13.43 | – | 94.30 | 300 | 0.854 | 1809 | 7/9 |
| NN+RL rl4L | 9 | 68.2 | 80.1 | 13.43 | 13.43 | – | 94.30 | 300 | 0.839 | 2391 | 6/9 |
| NN+RL rl4n29f | 9 | 70.3 | 75.5 | 13.83 | 13.83 | – | 94.13 | 300 | 0.860 | 2137 | 6/9 |

### smallval: protocol top300, reference ccsdt

| candidate | n | LUCJ % CCSD corr | LUCJ err (mHa) | QSCI err seed 0 | QSCI err 5-seed | seed sd | QSCI % corr | dim/spin | p_HF | configs in 10^5 | QSCI better than label |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 9 | 16.6 | 195.4 | 47.47 | 47.47 | – | 79.53 | 167 | 0.968 | 173 | 0/9 |
| compressed DF (Lin) | 9 | -401.3 | 1125.6 | 16.21 | 16.21 | – | 93.06 | 300 | 0.320 | 10433 | 0/9 |
| optimize=True label | 9 | 62.0 | 92.8 | 11.98 | 11.98 | – | 94.86 | 300 | 0.842 | 2097 | — |
| pretrained x4 | 9 | 59.4 | 97.6 | 11.52 | 11.52 | – | 95.05 | 300 | 0.854 | 1809 | 7/9 |
| NN+RL rl4L | 9 | 68.2 | 78.2 | 11.52 | 11.52 | – | 95.05 | 300 | 0.839 | 2391 | 6/9 |
| NN+RL rl4n29f | 9 | 70.3 | 73.6 | 11.92 | 11.92 | – | 94.88 | 300 | 0.860 | 2137 | 6/9 |

### smallval: protocol top300, candidate vs label (QSCI: 5-seed means)

| candidate vs optimize=True label | n | mean dE_LUCJ (mHa) | LUCJ lower | mean dE_QSCI (mHa) | QSCI lower | QSCI lower by > 2 s.e. | QSCI higher by > 2 s.e. |
|:--|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 9 | +102.6 | 0/9 | +35.49 | 0/9 | 0/9 | 9/9 |
| compressed DF (Lin) | 9 | +1032.7 | 0/9 | +4.23 | 0/9 | 0/9 | 9/9 |
| pretrained x4 | 9 | +4.8 | 5/9 | -0.46 | 7/9 | 7/9 | 2/9 |
| NN+RL rl4L | 9 | -14.6 | 8/9 | -0.46 | 6/9 | 6/9 | 3/9 |
| NN+RL rl4n29f | 9 | -19.2 | 9/9 | -0.06 | 6/9 | 6/9 | 3/9 |

* ranking LUCJ vs QSCI (all 6): rho per molecule mean 0.66 (median 0.77, min -0.09, max 0.89), pooled 0.14; concordant pairs 102/135 (resolved beyond 2 s.e.: 102/135); per molecule [0.89, 0.77, 0.54, -0.09, 0.89, 0.77, 0.77, 0.89, 0.54]
  pooled Spearman of QSCI energy with dim/spin -0.732, p_HF 0.448, entropy -0.461, LUCJ energy 0.14 (n=54)
* ranking LUCJ vs QSCI (without cdf_lin): rho per molecule mean 0.61 (median 0.70, min 0.10, max 0.90), pooled 0.45; concordant pairs 68/90 (resolved beyond 2 s.e.: 68/90); per molecule [0.9, 0.7, 0.3, 0.1, 0.9, 0.7, 0.7, 0.9, 0.3]
  pooled Spearman of QSCI energy with dim/spin -0.837, p_HF 0.529, entropy -0.5, LUCJ energy 0.452 (n=45)
* ranking LUCJ vs QSCI (label + network family): rho per molecule mean 0.22 (median 0.40, min -0.80, max 0.80), pooled -0.15; concordant pairs 32/54 (resolved beyond 2 s.e.: 32/54); per molecule [0.8, 0.4, -0.4, -0.8, 0.8, 0.4, 0.4, 0.8, -0.4]
  pooled Spearman of QSCI energy with dim/spin nan, p_HF 0.488, entropy -0.426, LUCJ energy -0.147 (n=36)

### smallval: per molecule, QSCI error in mHa (protocol lin, reference fci, mean of the draws) and LUCJ % CCSD corr in parentheses

| molecule | truncated CCSD | compressed DF (Lin) | optimize=True label | pretrained x4 | NN+RL rl4L | NN+RL rl4n29f | best QSCI | best LUCJ |
|:--|--:|--:|--:|--:|--:|--:|:--|:--|
| C2H3N_rxn2858_P | 65.0 (13) | 10.4 (-976) | 8.4 (62) | 8.4 (63) | 6.3 (66) | 7.3 (68) | NN+RL rl4L | NN+RL rl4n29f |
| C2H3N_rxn2858_TS | 77.4 (25) | 5.3 (-96) | 8.6 (57) | 10.1 (62) | 7.4 (69) | 8.7 (70) | compressed DF (Lin) | NN+RL rl4n29f |
| C2H4O_rxn0724_P | 73.5 (14) | 5.4 (-396) | 9.5 (67) | 9.7 (42) | 8.0 (66) | 10.4 (67) | compressed DF (Lin) | NN+RL rl4n29f |
| C2H4O_rxn0724_TS | 79.7 (12) | 8.8 (-326) | 11.3 (53) | 13.6 (42) | 13.1 (66) | 14.5 (66) | compressed DF (Lin) | NN+RL rl4L |
| C2H4O_rxn2507_P | 98.7 (18) | 4.4 (-517) | 8.1 (69) | 10.1 (59) | 7.7 (70) | 8.5 (73) | compressed DF (Lin) | NN+RL rl4n29f |
| C2H4O_rxn2507_TS | 60.3 (15) | 12.6 (-564) | 11.9 (49) | 14.2 (65) | 11.1 (68) | 13.2 (69) | NN+RL rl4L | NN+RL rl4n29f |
| C3H4_rxn2391_R | 97.6 (24) | 5.9 (-61) | 11.2 (70) | 11.5 (71) | 7.0 (75) | 8.7 (78) | compressed DF (Lin) | NN+RL rl4n29f |
| C3H4_rxn2391_P | 109.0 (16) | 6.9 (-177) | 9.0 (69) | 10.2 (67) | 7.2 (69) | 7.8 (72) | compressed DF (Lin) | NN+RL rl4n29f |
| C3H4_rxn2391_TS | 103.9 (13) | 12.4 (-497) | 12.7 (61) | 14.9 (64) | 9.8 (65) | 10.8 (69) | NN+RL rl4L | NN+RL rl4n29f |

### smallval: per molecule and candidate (protocol lin, reference fci, mHa)

| molecule | candidate | LUCJ % CCSD corr | LUCJ err | QSCI err seed 0: mean (min / max) | distinct batches | 5-seed mean ± sd | QSCI % corr | dim/spin | p_HF |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| C2H3N_rxn2858_P | truncated CCSD | 13.1 | 199.6 | 65.18 (65.18 / 65.18) | 1 | 64.96 ± 5.86 (5) | 71.60 | 82 | 0.985 |
| C2H3N_rxn2858_P | compressed DF (Lin) | -976.4 | 2402.0 | 10.21 (8.73 / 12.91) | 10 | 10.38 ± 0.29 (5) | 95.46 | 839 | 0.054 |
| C2H3N_rxn2858_P | optimize=True label | 61.8 | 91.1 | 8.12 (8.12 / 8.12) | 1 | 8.43 ± 0.58 (5) | 96.32 | 448 | 0.868 |
| C2H3N_rxn2858_P | pretrained x4 | 62.5 | 89.5 | 8.19 (8.19 / 8.19) | 1 | 8.40 ± 0.39 (5) | 96.33 | 386 | 0.871 |
| C2H3N_rxn2858_P | NN+RL rl4L | 65.5 | 82.9 | 6.36 (6.36 / 6.36) | 1 | 6.32 ± 0.27 (5) | 97.24 | 520 | 0.853 |
| C2H3N_rxn2858_P | NN+RL rl4n29f | 68.1 | 77.2 | 7.27 (7.27 / 7.27) | 1 | 7.25 ± 0.18 (5) | 96.83 | 465 | 0.873 |
| C2H3N_rxn2858_TS | truncated CCSD | 25.4 | 209.8 | 83.12 (83.12 / 83.12) | 1 | 77.43 ± 6.46 (5) | 72.02 | 58 | 0.958 |
| C2H3N_rxn2858_TS | compressed DF (Lin) | -96.1 | 530.2 | 5.40 (4.97 / 5.84) | 10 | 5.27 ± 0.18 (5) | 98.09 | 828 | 0.453 |
| C2H3N_rxn2858_TS | optimize=True label | 56.7 | 127.2 | 8.97 (8.97 / 8.97) | 1 | 8.61 ± 0.33 (5) | 96.89 | 489 | 0.794 |
| C2H3N_rxn2858_TS | pretrained x4 | 62.4 | 112.1 | 10.10 (10.10 / 10.10) | 1 | 10.11 ± 0.47 (5) | 96.35 | 407 | 0.794 |
| C2H3N_rxn2858_TS | NN+RL rl4L | 69.1 | 94.5 | 7.02 (7.02 / 7.02) | 1 | 7.39 ± 0.27 (5) | 97.33 | 551 | 0.792 |
| C2H3N_rxn2858_TS | NN+RL rl4n29f | 70.3 | 91.5 | 8.44 (8.44 / 8.44) | 1 | 8.72 ± 0.58 (5) | 96.85 | 508 | 0.819 |
| C2H4O_rxn0724_P | truncated CCSD | 14.1 | 190.2 | 69.79 (69.79 / 69.79) | 1 | 73.53 ± 9.55 (5) | 66.58 | 90 | 0.975 |
| C2H4O_rxn0724_P | compressed DF (Lin) | -395.7 | 1058.5 | 5.47 (5.13 / 5.83) | 10 | 5.36 ± 0.14 (5) | 97.56 | 844 | 0.370 |
| C2H4O_rxn0724_P | optimize=True label | 67.3 | 77.3 | 8.83 (8.83 / 8.83) | 1 | 9.48 ± 0.39 (5) | 95.69 | 468 | 0.879 |
| C2H4O_rxn0724_P | pretrained x4 | 41.8 | 131.5 | 9.39 (9.39 / 9.39) | 1 | 9.71 ± 0.29 (5) | 95.58 | 454 | 0.855 |
| C2H4O_rxn0724_P | NN+RL rl4L | 66.0 | 80.2 | 7.79 (7.79 / 7.79) | 1 | 8.04 ± 0.27 (5) | 96.35 | 570 | 0.862 |
| C2H4O_rxn0724_P | NN+RL rl4n29f | 67.4 | 77.2 | 10.95 (10.95 / 10.95) | 1 | 10.39 ± 0.59 (5) | 95.28 | 458 | 0.886 |
| C2H4O_rxn0724_TS | truncated CCSD | 11.9 | 202.9 | 80.75 (80.75 / 80.75) | 1 | 79.73 ± 3.09 (5) | 65.22 | 69 | 0.970 |
| C2H4O_rxn0724_TS | compressed DF (Lin) | -326.0 | 951.5 | 8.57 (7.91 / 9.10) | 10 | 8.75 ± 0.25 (5) | 96.18 | 890 | 0.289 |
| C2H4O_rxn0724_TS | optimize=True label | 53.3 | 111.2 | 11.77 (11.77 / 11.77) | 1 | 11.32 ± 0.43 (5) | 95.06 | 595 | 0.826 |
| C2H4O_rxn0724_TS | pretrained x4 | 42.2 | 135.9 | 14.20 (14.20 / 14.20) | 1 | 13.59 ± 0.37 (5) | 94.07 | 491 | 0.860 |
| C2H4O_rxn0724_TS | NN+RL rl4L | 66.3 | 82.4 | 12.96 (12.96 / 12.96) | 1 | 13.14 ± 0.30 (5) | 94.27 | 552 | 0.877 |
| C2H4O_rxn0724_TS | NN+RL rl4n29f | 66.0 | 83.0 | 13.86 (13.86 / 13.86) | 1 | 14.45 ± 0.62 (5) | 93.70 | 515 | 0.892 |
| C2H4O_rxn2507_P | truncated CCSD | 18.2 | 165.6 | 81.93 (81.93 / 81.93) | 1 | 98.73 ± 13.00 (5) | 51.05 | 67 | 0.986 |
| C2H4O_rxn2507_P | compressed DF (Lin) | -517.1 | 1224.7 | 4.48 (4.17 / 5.13) | 10 | 4.38 ± 0.10 (5) | 97.83 | 973 | 0.320 |
| C2H4O_rxn2507_P | optimize=True label | 69.4 | 64.4 | 7.61 (7.61 / 7.61) | 1 | 8.12 ± 0.38 (5) | 95.97 | 421 | 0.885 |
| C2H4O_rxn2507_P | pretrained x4 | 59.3 | 84.3 | 9.86 (9.86 / 9.86) | 1 | 10.05 ± 0.33 (5) | 95.02 | 358 | 0.897 |
| C2H4O_rxn2507_P | NN+RL rl4L | 70.0 | 63.3 | 7.84 (7.84 / 7.84) | 1 | 7.66 ± 0.21 (5) | 96.20 | 486 | 0.876 |
| C2H4O_rxn2507_P | NN+RL rl4n29f | 72.6 | 58.0 | 8.84 (8.84 / 8.84) | 1 | 8.52 ± 0.30 (5) | 95.77 | 423 | 0.896 |
| C2H4O_rxn2507_TS | truncated CCSD | 14.5 | 205.9 | 55.94 (55.94 / 55.94) | 1 | 60.32 ± 11.55 (5) | 74.75 | 168 | 0.972 |
| C2H4O_rxn2507_TS | compressed DF (Lin) | -564.0 | 1522.3 | 12.86 (11.54 / 13.94) | 10 | 12.63 ± 0.30 (5) | 94.71 | 1052 | 0.131 |
| C2H4O_rxn2507_TS | optimize=True label | 49.5 | 126.3 | 11.24 (11.24 / 11.24) | 1 | 11.91 ± 0.68 (5) | 95.02 | 605 | 0.845 |
| C2H4O_rxn2507_TS | pretrained x4 | 64.6 | 91.9 | 14.15 (14.15 / 14.15) | 1 | 14.18 ± 0.11 (5) | 94.06 | 509 | 0.872 |
| C2H4O_rxn2507_TS | NN+RL rl4L | 67.5 | 85.3 | 10.30 (10.30 / 10.30) | 1 | 11.05 ± 0.44 (5) | 95.37 | 656 | 0.866 |
| C2H4O_rxn2507_TS | NN+RL rl4n29f | 69.3 | 81.1 | 13.07 (13.07 / 13.07) | 1 | 13.17 ± 0.62 (5) | 94.49 | 557 | 0.882 |
| C3H4_rxn2391_R | truncated CCSD | 24.0 | 180.1 | 88.71 (88.71 / 88.71) | 1 | 97.57 ± 11.55 (5) | 58.66 | 72 | 0.985 |
| C3H4_rxn2391_R | compressed DF (Lin) | -61.5 | 379.0 | 5.87 (5.60 / 6.20) | 10 | 5.91 ± 0.39 (5) | 97.50 | 921 | 0.672 |
| C3H4_rxn2391_R | optimize=True label | 70.1 | 73.0 | 10.80 (10.80 / 10.80) | 1 | 11.23 ± 0.32 (5) | 95.24 | 432 | 0.885 |
| C3H4_rxn2391_R | pretrained x4 | 71.2 | 70.3 | 11.87 (11.87 / 11.87) | 1 | 11.52 ± 0.60 (5) | 95.12 | 411 | 0.897 |
| C3H4_rxn2391_R | NN+RL rl4L | 75.5 | 60.4 | 6.87 (6.87 / 6.87) | 1 | 6.99 ± 0.15 (5) | 97.04 | 648 | 0.857 |
| C3H4_rxn2391_R | NN+RL rl4n29f | 77.8 | 54.9 | 8.21 (8.21 / 8.21) | 1 | 8.68 ± 0.48 (5) | 96.32 | 535 | 0.881 |
| C3H4_rxn2391_P | truncated CCSD | 15.6 | 201.1 | 112.10 (112.10 / 112.10) | 1 | 109.03 ± 7.39 (5) | 54.11 | 54 | 0.985 |
| C3H4_rxn2391_P | compressed DF (Lin) | -177.1 | 652.8 | 6.75 (6.50 / 7.00) | 10 | 6.92 ± 0.41 (5) | 97.09 | 894 | 0.467 |
| C3H4_rxn2391_P | optimize=True label | 68.8 | 76.3 | 8.59 (8.59 / 8.59) | 1 | 9.04 ± 0.76 (5) | 96.20 | 545 | 0.851 |
| C3H4_rxn2391_P | pretrained x4 | 66.7 | 81.3 | 10.09 (10.09 / 10.09) | 1 | 10.23 ± 0.38 (5) | 95.69 | 426 | 0.868 |
| C3H4_rxn2391_P | NN+RL rl4L | 69.2 | 75.5 | 7.08 (7.08 / 7.08) | 1 | 7.18 ± 0.13 (5) | 96.98 | 641 | 0.832 |
| C3H4_rxn2391_P | NN+RL rl4n29f | 71.9 | 68.9 | 7.99 (7.99 / 7.99) | 1 | 7.79 ± 0.41 (5) | 96.72 | 565 | 0.849 |
| C3H4_rxn2391_TS | truncated CCSD | 12.5 | 221.0 | 103.76 (103.76 / 103.76) | 1 | 103.85 ± 3.19 (5) | 58.57 | 90 | 0.900 |
| C3H4_rxn2391_TS | compressed DF (Lin) | -497.4 | 1426.4 | 12.41 (11.49 / 13.57) | 10 | 12.45 ± 0.15 (5) | 95.03 | 1118 | 0.129 |
| C3H4_rxn2391_TS | optimize=True label | 61.2 | 106.0 | 12.32 (12.32 / 12.32) | 1 | 12.74 ± 0.28 (5) | 94.92 | 679 | 0.743 |
| C3H4_rxn2391_TS | pretrained x4 | 64.1 | 99.0 | 15.08 (15.08 / 15.08) | 1 | 14.91 ± 0.79 (5) | 94.05 | 593 | 0.774 |
| C3H4_rxn2391_TS | NN+RL rl4L | 65.2 | 96.6 | 9.44 (9.44 / 9.44) | 1 | 9.76 ± 0.26 (5) | 96.11 | 875 | 0.739 |
| C3H4_rxn2391_TS | NN+RL rl4n29f | 68.9 | 87.8 | 10.70 (10.70 / 10.70) | 1 | 10.80 ± 0.27 (5) | 95.69 | 816 | 0.759 |

### smallval: costs

| stage | median s | n |
|:--|--:|--:|
| t_energy (per state / per draw) | 2.2 | 270 |
| t_sample (per state / per draw) | 3.5 | 270 |
| t_ci (per state / per draw) | 19.0 | 270 |
| diag lin, norb 15 (threads [2, 4]) | 19.4 (max 37) | 50 |
| diag lin, norb 15, cdf_lin (threads [4]) | 30.0 (max 41) | 100 |
| diag lin, norb 16 (threads [4]) | 33.1 (max 82) | 175 |
| diag lin, norb 16, cdf_lin (threads [4]) | 55.9 (max 88) | 350 |
| diag n2631g, norb 15 (threads [4]) | 54.4 (max 105) | 30 |
| diag n2631g, norb 15, cdf_lin (threads [4]) | 196.7 (max 248) | 2 |
| diag n2631g, norb 16 (threads [4]) | 100.2 (max 310) | 105 |
| diag n2631g, norb 16, cdf_lin (threads [4]) | 414.4 (max 1009) | 7 |
| diag top100, norb 15 (threads [2]) | 2.4 (max 3) | 10 |
| diag top100, norb 15, cdf_lin (threads [2]) | 2.2 (max 3) | 2 |
| diag top100, norb 16 (threads [2]) | 3.0 (max 5) | 35 |
| diag top100, norb 16, cdf_lin (threads [2]) | 3.2 (max 5) | 7 |
| diag top300, norb 15 (threads [2]) | 14.7 (max 18) | 10 |
| diag top300, norb 15, cdf_lin (threads [2]) | 13.3 (max 15) | 2 |
| diag top300, norb 16 (threads [2]) | 22.0 (max 33) | 35 |
| diag top300, norb 16, cdf_lin (threads [2]) | 24.6 (max 27) | 7 |
