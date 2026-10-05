
### n17: protocol lin, reference ccsdt

| candidate | n | LUCJ % CCSD corr | LUCJ err (mHa) | QSCI err seed 0 | QSCI err 5-seed | seed sd | QSCI % corr | dim/spin | p_HF | configs in 10^5 | QSCI better than label |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 12 | 15.5 | 249.1 | 79.88 | 83.59 | 5.29 | 71.31 | 86 | 0.973 | 177 | 0/12 |
| compressed DF (Lin) | 12 | -145.1 | 704.7 | 8.42 | 8.27 | 0.31 | 97.16 | 1021 | 0.429 | 8299 | 12/12 |
| optimize=True label | 12 | 61.1 | 119.0 | 10.39 | 10.42 | 0.51 | 96.43 | 594 | 0.829 | 2170 | — |
| pretrained x4 | 12 | 51.9 | 145.1 | 10.37 | 10.42 | 0.44 | 96.42 | 560 | 0.823 | 2029 | 7/12 |
| NN+RL rl4L | 12 | 69.9 | 94.1 | 9.08 | 8.89 | 0.38 | 96.94 | 654 | 0.842 | 2258 | 12/12 |
| NN+RL rl4n29f | 12 | 68.4 | 98.1 | 11.66 | 11.69 | 0.60 | 95.99 | 523 | 0.865 | 1859 | 0/12 |

QSCI error ratio truncated / candidate (5-seed means, reference ccsdt):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| compressed DF (Lin) | 12 | 11.0 | 10.5 | 6.2 | 20.4 |
| optimize=True label | 12 | 8.4 | 8.2 | 6.1 | 14.5 |
| pretrained x4 | 12 | 8.3 | 8.1 | 6.1 | 12.2 |
| NN+RL rl4L | 12 | 9.8 | 9.6 | 6.7 | 15.9 |
| NN+RL rl4n29f | 12 | 7.5 | 7.3 | 5.2 | 12.9 |

QSCI error ratio optimize=True label / candidate (reference ccsdt; > 1 = candidate better; molecules with an error <= 0 skipped):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| compressed DF (Lin) | 12 | 1.29 | 1.28 | 1.03 | 1.45 |
| pretrained x4 | 12 | 1.00 | 1.00 | 0.84 | 1.14 |
| NN+RL rl4L | 12 | 1.18 | 1.17 | 1.05 | 1.36 |
| NN+RL rl4n29f | 12 | 0.89 | 0.89 | 0.79 | 1.00 |

QSCI error ratio compressed DF (Lin) / candidate (reference ccsdt; > 1 = candidate better; molecules with an error <= 0 skipped):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| optimize=True label | 12 | 0.78 | 0.78 | 0.69 | 0.97 |
| pretrained x4 | 12 | 0.78 | 0.78 | 0.60 | 1.07 |
| NN+RL rl4L | 12 | 0.92 | 0.91 | 0.77 | 1.32 |
| NN+RL rl4n29f | 12 | 0.70 | 0.70 | 0.60 | 0.97 |

### n17: protocol lin, candidate vs label (QSCI: 5-seed means)

| candidate vs optimize=True label | n | mean dE_LUCJ (mHa) | LUCJ lower | mean dE_QSCI (mHa) | QSCI lower | QSCI lower by > 2 s.e. | QSCI higher by > 2 s.e. |
|:--|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 12 | +130.0 | 0/12 | +73.17 | 0/12 | 0/12 | 12/12 |
| compressed DF (Lin) | 12 | +585.7 | 0/12 | -2.14 | 12/12 | 11/12 | 0/12 |
| pretrained x4 | 12 | +26.1 | 1/12 | +0.01 | 7/12 | 2/12 | 2/12 |
| NN+RL rl4L | 12 | -25.0 | 10/12 | -1.52 | 12/12 | 10/12 | 0/12 |
| NN+RL rl4n29f | 12 | -20.9 | 9/12 | +1.27 | 0/12 | 0/12 | 11/12 |

* ranking LUCJ vs QSCI (all 6): rho per molecule mean -0.02 (median -0.09, min -0.26, max 0.31), pooled 0.06; concordant pairs 91/180 (resolved beyond 2 s.e.: 85/162); per molecule [-0.03, -0.09, 0.31, 0.09, -0.09, -0.09, 0.26, -0.09, -0.26, 0.09, -0.09, -0.26]
  pooled Spearman of QSCI energy with dim/spin -0.781, p_HF 0.571, entropy -0.604, LUCJ energy 0.061 (n=72)
* ranking LUCJ vs QSCI (without cdf_lin): rho per molecule mean 0.62 (median 0.60, min 0.30, max 0.90), pooled 0.50; concordant pairs 88/120 (resolved beyond 2 s.e.: 83/106); per molecule [0.7, 0.6, 0.6, 0.9, 0.6, 0.6, 0.8, 0.6, 0.3, 0.9, 0.6, 0.3]
  pooled Spearman of QSCI energy with dim/spin -0.743, p_HF 0.534, entropy -0.537, LUCJ energy 0.502 (n=60)
* ranking LUCJ vs QSCI (label + network family): rho per molecule mean 0.25 (median 0.20, min -0.40, max 0.80), pooled 0.12; concordant pairs 40/72 (resolved beyond 2 s.e.: 35/58); per molecule [0.4, 0.2, 0.2, 0.8, 0.2, 0.2, 0.6, 0.2, -0.4, 0.8, 0.2, -0.4]
  pooled Spearman of QSCI energy with dim/spin -0.893, p_HF 0.379, entropy -0.414, LUCJ energy 0.121 (n=48)

### n17: protocol n2631g, reference ccsdt

| candidate | n | LUCJ % CCSD corr | LUCJ err (mHa) | QSCI err seed 0 | QSCI err 5-seed | seed sd | QSCI % corr | dim/spin | p_HF | configs in 10^5 | QSCI better than label |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 12 | 15.5 | 249.1 | 52.10 | 52.10 | – | 82.15 | 183 | 0.973 | 177 | 0/12 |
| compressed DF (Lin) | 12 | -145.1 | 704.7 | -1.62 | -1.62 | – | 100.54 | 3787 | 0.429 | 8299 | 12/12 |
| optimize=True label | 12 | 61.1 | 119.0 | 0.94 | 0.94 | – | 99.67 | 1478 | 0.829 | 2170 | — |
| pretrained x4 | 12 | 51.9 | 145.1 | 1.16 | 1.16 | – | 99.59 | 1388 | 0.823 | 2029 | 3/12 |
| NN+RL rl4L | 12 | 69.9 | 94.1 | 0.29 | 0.29 | – | 99.89 | 1666 | 0.842 | 2258 | 11/12 |
| NN+RL rl4n29f | 12 | 68.4 | 98.1 | 1.56 | 1.56 | – | 99.45 | 1324 | 0.865 | 1859 | 1/12 |

QSCI error ratio truncated / candidate (5-seed means, reference ccsdt):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| optimize=True label | 10 | 50.7 | 47.1 | 22.4 | 65.8 |
| pretrained x4 | 10 | 36.9 | 35.7 | 18.8 | 43.8 |
| NN+RL rl4L | 6 | 314.6 | 92.1 | 26.8 | 1596.0 |
| NN+RL rl4n29f | 11 | 40.0 | 33.8 | 15.4 | 124.7 |

QSCI error ratio optimize=True label / candidate (reference ccsdt; > 1 = candidate better; molecules with an error <= 0 skipped):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|
| pretrained x4 | 10 | 0.78 | 0.76 | 0.59 | 1.29 |
| NN+RL rl4L | 6 | 9.15 | 2.44 | 0.79 | 47.65 |
| NN+RL rl4n29f | 10 | 0.65 | 0.63 | 0.37 | 1.11 |

QSCI error ratio compressed DF (Lin) / candidate (reference ccsdt; > 1 = candidate better; molecules with an error <= 0 skipped):

| candidate | n | mean | geometric mean | min | max |
|:--|--:|--:|--:|--:|--:|

### n17: protocol n2631g, candidate vs label (QSCI: 5-seed means)

| candidate vs optimize=True label | n | mean dE_LUCJ (mHa) | LUCJ lower | mean dE_QSCI (mHa) | QSCI lower | QSCI lower by > 2 s.e. | QSCI higher by > 2 s.e. |
|:--|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 12 | +130.0 | 0/12 | +51.16 | 0/12 | 0/12 | 12/12 |
| compressed DF (Lin) | 12 | +585.7 | 0/12 | -2.55 | 12/12 | 12/12 | 0/12 |
| pretrained x4 | 12 | +26.1 | 1/12 | +0.23 | 3/12 | 3/12 | 9/12 |
| NN+RL rl4L | 12 | -25.0 | 10/12 | -0.65 | 11/12 | 11/12 | 1/12 |
| NN+RL rl4n29f | 12 | -20.9 | 9/12 | +0.62 | 1/12 | 1/12 | 11/12 |

* ranking LUCJ vs QSCI (all 6): rho per molecule mean -0.05 (median -0.03, min -0.26, max 0.09), pooled 0.05; concordant pairs 93/180 (resolved beyond 2 s.e.: 93/180); per molecule [-0.03, -0.03, -0.03, 0.03, -0.03, -0.09, 0.03, -0.03, -0.2, 0.09, -0.03, -0.26]
  pooled Spearman of QSCI energy with dim/spin -0.779, p_HF 0.652, entropy -0.661, LUCJ energy 0.05 (n=72)
* ranking LUCJ vs QSCI (without cdf_lin): rho per molecule mean 0.67 (median 0.70, min 0.30, max 0.90), pooled 0.54; concordant pairs 93/120 (resolved beyond 2 s.e.: 93/120); per molecule [0.7, 0.7, 0.7, 0.8, 0.7, 0.6, 0.8, 0.7, 0.4, 0.9, 0.7, 0.3]
  pooled Spearman of QSCI energy with dim/spin -0.716, p_HF 0.542, entropy -0.53, LUCJ energy 0.536 (n=60)
* ranking LUCJ vs QSCI (label + network family): rho per molecule mean 0.33 (median 0.40, min -0.40, max 0.80), pooled 0.18; concordant pairs 45/72 (resolved beyond 2 s.e.: 45/72); per molecule [0.4, 0.4, 0.4, 0.6, 0.4, 0.2, 0.6, 0.4, -0.2, 0.8, 0.4, -0.4]
  pooled Spearman of QSCI energy with dim/spin -0.959, p_HF 0.351, entropy -0.429, LUCJ energy 0.177 (n=48)

### n17: protocol top100, reference ccsdt

| candidate | n | LUCJ % CCSD corr | LUCJ err (mHa) | QSCI err seed 0 | QSCI err 5-seed | seed sd | QSCI % corr | dim/spin | p_HF | configs in 10^5 | QSCI better than label |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 12 | 15.5 | 249.1 | 67.85 | 67.85 | – | 76.76 | 98 | 0.973 | 177 | 0/12 |
| compressed DF (Lin) | 12 | -145.1 | 704.7 | 45.74 | 45.74 | – | 84.35 | 100 | 0.429 | 8299 | 0/12 |
| optimize=True label | 12 | 61.1 | 119.0 | 40.54 | 40.54 | – | 86.15 | 100 | 0.829 | 2170 | — |
| pretrained x4 | 12 | 51.9 | 145.1 | 39.80 | 39.80 | – | 86.42 | 100 | 0.823 | 2029 | 8/12 |
| NN+RL rl4L | 12 | 69.9 | 94.1 | 39.41 | 39.41 | – | 86.54 | 100 | 0.842 | 2258 | 10/12 |
| NN+RL rl4n29f | 12 | 68.4 | 98.1 | 39.96 | 39.96 | – | 86.35 | 100 | 0.865 | 1859 | 8/12 |

### n17: protocol top100, candidate vs label (QSCI: 5-seed means)

| candidate vs optimize=True label | n | mean dE_LUCJ (mHa) | LUCJ lower | mean dE_QSCI (mHa) | QSCI lower | QSCI lower by > 2 s.e. | QSCI higher by > 2 s.e. |
|:--|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 12 | +130.0 | 0/12 | +27.31 | 0/12 | 0/12 | 12/12 |
| compressed DF (Lin) | 12 | +585.7 | 0/12 | +5.20 | 0/12 | 0/12 | 12/12 |
| pretrained x4 | 12 | +26.1 | 1/12 | -0.74 | 8/12 | 8/12 | 4/12 |
| NN+RL rl4L | 12 | -25.0 | 10/12 | -1.13 | 10/12 | 10/12 | 2/12 |
| NN+RL rl4n29f | 12 | -20.9 | 9/12 | -0.57 | 8/12 | 8/12 | 4/12 |

* ranking LUCJ vs QSCI (all 6): rho per molecule mean 0.68 (median 0.77, min 0.43, max 0.77), pooled 0.69; concordant pairs 135/180 (resolved beyond 2 s.e.: 135/180); per molecule [0.43, 0.77, 0.6, 0.77, 0.77, 0.77, 0.43, 0.77, 0.66, 0.6, 0.77, 0.77]
  pooled Spearman of QSCI energy with dim/spin -0.229, p_HF 0.018, entropy -0.054, LUCJ energy 0.691 (n=72)
* ranking LUCJ vs QSCI (without cdf_lin): rho per molecule mean 0.53 (median 0.70, min 0.10, max 0.70), pooled 0.46; concordant pairs 87/120 (resolved beyond 2 s.e.: 87/120); per molecule [0.1, 0.7, 0.4, 0.7, 0.7, 0.7, 0.1, 0.7, 0.5, 0.4, 0.7, 0.7]
  pooled Spearman of QSCI energy with dim/spin -0.302, p_HF 0.385, entropy -0.348, LUCJ energy 0.463 (n=60)
* ranking LUCJ vs QSCI (label + network family): rho per molecule mean 0.07 (median 0.40, min -0.80, max 0.40), pooled 0.13; concordant pairs 39/72 (resolved beyond 2 s.e.: 39/72); per molecule [-0.8, 0.4, -0.2, 0.4, 0.4, 0.4, -0.8, 0.4, 0.0, -0.2, 0.4, 0.4]
  pooled Spearman of QSCI energy with dim/spin nan, p_HF -0.066, entropy 0.152, LUCJ energy 0.132 (n=48)

### n17: protocol top300, reference ccsdt

| candidate | n | LUCJ % CCSD corr | LUCJ err (mHa) | QSCI err seed 0 | QSCI err 5-seed | seed sd | QSCI % corr | dim/spin | p_HF | configs in 10^5 | QSCI better than label |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 12 | 15.5 | 249.1 | 52.10 | 52.10 | – | 82.15 | 183 | 0.973 | 177 | 0/12 |
| compressed DF (Lin) | 12 | -145.1 | 704.7 | 22.18 | 22.18 | – | 92.41 | 300 | 0.429 | 8299 | 1/12 |
| optimize=True label | 12 | 61.1 | 119.0 | 16.86 | 16.86 | – | 94.24 | 300 | 0.829 | 2170 | — |
| pretrained x4 | 12 | 51.9 | 145.1 | 15.59 | 15.59 | – | 94.67 | 300 | 0.823 | 2029 | 11/12 |
| NN+RL rl4L | 12 | 69.9 | 94.1 | 15.74 | 15.74 | – | 94.62 | 300 | 0.842 | 2258 | 11/12 |
| NN+RL rl4n29f | 12 | 68.4 | 98.1 | 16.15 | 16.15 | – | 94.48 | 300 | 0.865 | 1859 | 10/12 |

### n17: protocol top300, candidate vs label (QSCI: 5-seed means)

| candidate vs optimize=True label | n | mean dE_LUCJ (mHa) | LUCJ lower | mean dE_QSCI (mHa) | QSCI lower | QSCI lower by > 2 s.e. | QSCI higher by > 2 s.e. |
|:--|--:|--:|--:|--:|--:|--:|--:|
| truncated CCSD | 12 | +130.0 | 0/12 | +35.24 | 0/12 | 0/12 | 12/12 |
| compressed DF (Lin) | 12 | +585.7 | 0/12 | +5.32 | 1/12 | 1/12 | 11/12 |
| pretrained x4 | 12 | +26.1 | 1/12 | -1.26 | 11/12 | 11/12 | 1/12 |
| NN+RL rl4L | 12 | -25.0 | 10/12 | -1.12 | 11/12 | 11/12 | 1/12 |
| NN+RL rl4n29f | 12 | -20.9 | 9/12 | -0.71 | 10/12 | 10/12 | 2/12 |

* ranking LUCJ vs QSCI (all 6): rho per molecule mean 0.59 (median 0.60, min -0.14, max 0.94), pooled 0.74; concordant pairs 130/180 (resolved beyond 2 s.e.: 130/180); per molecule [0.54, 0.6, 0.54, 0.94, 0.6, 0.77, 0.54, 0.6, 0.77, 0.71, 0.6, -0.14]
  pooled Spearman of QSCI energy with dim/spin -0.717, p_HF 0.103, entropy -0.067, LUCJ energy 0.74 (n=72)
* ranking LUCJ vs QSCI (without cdf_lin): rho per molecule mean 0.48 (median 0.40, min 0.30, max 1.00), pooled 0.44; concordant pairs 85/120 (resolved beyond 2 s.e.: 85/120); per molecule [0.3, 0.4, 0.3, 1.0, 0.4, 0.7, 0.3, 0.4, 0.7, 0.6, 0.4, 0.3]
  pooled Spearman of QSCI energy with dim/spin -0.911, p_HF 0.474, entropy -0.41, LUCJ energy 0.438 (n=60)
* ranking LUCJ vs QSCI (label + network family): rho per molecule mean -0.03 (median -0.20, min -0.40, max 1.00), pooled -0.14; concordant pairs 37/72 (resolved beyond 2 s.e.: 37/72); per molecule [-0.4, -0.2, -0.4, 1.0, -0.2, 0.4, -0.4, -0.2, 0.4, 0.2, -0.2, -0.4]
  pooled Spearman of QSCI energy with dim/spin nan, p_HF 0.186, entropy -0.019, LUCJ energy -0.138 (n=48)

### n17: per molecule, QSCI error in mHa (protocol lin, reference ccsdt, mean of the draws) and LUCJ % CCSD corr in parentheses

| molecule | truncated CCSD | compressed DF (Lin) | optimize=True label | pretrained x4 | NN+RL rl4L | NN+RL rl4n29f | best QSCI | best LUCJ |
|:--|--:|--:|--:|--:|--:|--:|:--|:--|
| C2HNO_rxn4439_P | 89.9 (20) | 9.8 (-261) | 10.6 (55) | 12.0 (50) | 10.0 (70) | 13.3 (69) | compressed DF (Lin) | NN+RL rl4L |
| C2HNO_rxn4439_R | 81.4 (17) | 6.9 (-122) | 9.5 (62) | 9.4 (56) | 7.9 (73) | 10.3 (71) | compressed DF (Lin) | NN+RL rl4L |
| C2HNO_rxn4439_TS | 80.7 (9) | 12.9 (-178) | 13.3 (57) | 12.1 (57) | 9.8 (69) | 13.4 (68) | NN+RL rl4L | NN+RL rl4L |
| C2HNO_rxn4440_P | 84.5 (18) | 4.1 (-316) | 5.8 (72) | 6.9 (53) | 5.3 (69) | 6.6 (66) | compressed DF (Lin) | optimize=True label |
| C2HNO_rxn4440_R | 81.4 (17) | 6.9 (-122) | 9.5 (62) | 9.4 (56) | 7.9 (73) | 10.3 (71) | compressed DF (Lin) | NN+RL rl4L |
| C2HNO_rxn4440_TS | 83.3 (15) | 6.5 (-88) | 9.4 (67) | 8.2 (46) | 7.2 (69) | 10.8 (68) | compressed DF (Lin) | NN+RL rl4L |
| C2HNO_rxn4441_P | 89.8 (13) | 10.7 (-158) | 11.7 (66) | 12.2 (51) | 10.3 (66) | 13.1 (64) | NN+RL rl4L | optimize=True label |
| C2HNO_rxn4441_R | 81.4 (17) | 6.9 (-122) | 9.5 (62) | 9.4 (56) | 7.9 (73) | 10.3 (71) | compressed DF (Lin) | NN+RL rl4L |
| C2HNO_rxn4441_TS | 78.3 (10) | 10.2 (-80) | 12.6 (46) | 12.9 (47) | 11.3 (67) | 15.0 (68) | compressed DF (Lin) | NN+RL rl4n29f |
| C2HNO_rxn4442_P | 95.6 (20) | 8.6 (-129) | 11.3 (70) | 11.6 (50) | 9.9 (71) | 12.8 (68) | compressed DF (Lin) | NN+RL rl4L |
| C2HNO_rxn4442_R | 81.4 (17) | 6.9 (-122) | 9.5 (62) | 9.4 (56) | 7.9 (73) | 10.3 (71) | compressed DF (Lin) | NN+RL rl4L |
| C2HNO_rxn4442_TS | 75.5 (13) | 8.7 (-43) | 12.1 (54) | 11.7 (43) | 11.3 (64) | 14.1 (64) | compressed DF (Lin) | NN+RL rl4n29f |

### n17: per molecule and candidate (protocol lin, reference ccsdt, mHa)

| molecule | candidate | LUCJ % CCSD corr | LUCJ err | QSCI err seed 0: mean (min / max) | distinct batches | 5-seed mean ± sd | QSCI % corr | dim/spin | p_HF |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| C2HNO_rxn4439_P | truncated CCSD | 19.6 | 216.4 | 94.83 (94.83 / 94.83) | 1 | 89.89 ± 8.47 (5) | 66.39 | 95 | 0.969 |
| C2HNO_rxn4439_P | compressed DF (Lin) | -261.1 | 947.7 | 10.47 (9.86 / 11.33) | 10 | 9.76 ± 0.44 (5) | 96.35 | 1000 | 0.322 |
| C2HNO_rxn4439_P | optimize=True label | 54.7 | 125.0 | 10.16 (10.16 / 10.16) | 1 | 10.56 ± 0.30 (5) | 96.05 | 602 | 0.830 |
| C2HNO_rxn4439_P | pretrained x4 | 50.2 | 136.6 | 12.74 (12.74 / 12.74) | 1 | 11.97 ± 0.64 (5) | 95.53 | 562 | 0.839 |
| C2HNO_rxn4439_P | NN+RL rl4L | 70.4 | 83.9 | 11.02 (11.02 / 11.02) | 1 | 10.04 ± 0.65 (5) | 96.25 | 655 | 0.861 |
| C2HNO_rxn4439_P | NN+RL rl4n29f | 68.7 | 88.3 | 13.82 (13.82 / 13.82) | 1 | 13.30 ± 0.91 (5) | 95.03 | 536 | 0.878 |
| C2HNO_rxn4439_R | truncated CCSD | 17.1 | 243.1 | 75.32 (75.32 / 75.32) | 1 | 81.36 ± 3.47 (5) | 72.18 | 73 | 0.985 |
| C2HNO_rxn4439_R | compressed DF (Lin) | -122.2 | 644.1 | 7.38 (6.75 / 7.99) | 10 | 6.93 ± 0.26 (5) | 97.63 | 1076 | 0.504 |
| C2HNO_rxn4439_R | optimize=True label | 61.6 | 115.2 | 9.16 (9.16 / 9.16) | 1 | 9.53 ± 0.43 (5) | 96.74 | 554 | 0.833 |
| C2HNO_rxn4439_R | pretrained x4 | 56.4 | 130.2 | 8.94 (8.94 / 8.94) | 1 | 9.37 ± 0.30 (5) | 96.80 | 494 | 0.840 |
| C2HNO_rxn4439_R | NN+RL rl4L | 73.3 | 81.5 | 8.10 (8.10 / 8.10) | 1 | 7.89 ± 0.28 (5) | 97.30 | 608 | 0.853 |
| C2HNO_rxn4439_R | NN+RL rl4n29f | 71.4 | 86.8 | 10.39 (10.39 / 10.39) | 1 | 10.32 ± 0.36 (5) | 96.47 | 477 | 0.878 |
| C2HNO_rxn4439_TS | truncated CCSD | 9.3 | 292.2 | 80.31 (80.31 / 80.31) | 1 | 80.66 ± 4.24 (5) | 74.85 | 111 | 0.960 |
| C2HNO_rxn4439_TS | compressed DF (Lin) | -177.7 | 868.0 | 12.97 (11.04 / 13.91) | 10 | 12.95 ± 0.22 (5) | 95.96 | 1154 | 0.329 |
| C2HNO_rxn4439_TS | optimize=True label | 57.3 | 144.3 | 13.11 (13.11 / 13.11) | 1 | 13.33 ± 0.75 (5) | 95.85 | 730 | 0.803 |
| C2HNO_rxn4439_TS | pretrained x4 | 56.6 | 146.4 | 11.67 (11.67 / 11.67) | 1 | 12.06 ± 0.54 (5) | 96.24 | 705 | 0.775 |
| C2HNO_rxn4439_TS | NN+RL rl4L | 69.4 | 106.8 | 10.53 (10.53 / 10.53) | 1 | 9.83 ± 0.44 (5) | 96.93 | 826 | 0.805 |
| C2HNO_rxn4439_TS | NN+RL rl4n29f | 68.5 | 109.7 | 13.34 (13.34 / 13.34) | 1 | 13.37 ± 0.63 (5) | 95.83 | 669 | 0.829 |
| C2HNO_rxn4440_P | truncated CCSD | 17.5 | 251.1 | 79.79 (79.79 / 79.79) | 1 | 84.55 ± 7.84 (5) | 72.08 | 41 | 0.976 |
| C2HNO_rxn4440_P | compressed DF (Lin) | -316.3 | 1233.6 | 4.09 (3.03 / 5.07) | 10 | 4.14 ± 0.38 (5) | 98.63 | 934 | 0.310 |
| C2HNO_rxn4440_P | optimize=True label | 71.9 | 91.1 | 5.97 (5.97 / 5.97) | 1 | 5.81 ± 0.75 (5) | 98.08 | 443 | 0.851 |
| C2HNO_rxn4440_P | pretrained x4 | 53.3 | 146.0 | 6.35 (6.35 / 6.35) | 1 | 6.93 ± 0.49 (5) | 97.71 | 459 | 0.832 |
| C2HNO_rxn4440_P | NN+RL rl4L | 69.2 | 99.2 | 4.77 (4.77 / 4.77) | 1 | 5.31 ± 0.59 (5) | 98.25 | 523 | 0.850 |
| C2HNO_rxn4440_P | NN+RL rl4n29f | 66.2 | 107.8 | 6.67 (6.67 / 6.67) | 1 | 6.58 ± 0.40 (5) | 97.83 | 426 | 0.871 |
| C2HNO_rxn4440_R | truncated CCSD | 17.1 | 243.1 | 75.32 (75.32 / 75.32) | 1 | 81.36 ± 3.47 (5) | 72.18 | 73 | 0.985 |
| C2HNO_rxn4440_R | compressed DF (Lin) | -122.2 | 644.1 | 7.38 (6.75 / 7.99) | 10 | 6.93 ± 0.26 (5) | 97.63 | 1076 | 0.504 |
| C2HNO_rxn4440_R | optimize=True label | 61.6 | 115.2 | 9.16 (9.16 / 9.16) | 1 | 9.53 ± 0.43 (5) | 96.74 | 554 | 0.833 |
| C2HNO_rxn4440_R | pretrained x4 | 56.4 | 130.2 | 8.94 (8.94 / 8.94) | 1 | 9.37 ± 0.30 (5) | 96.80 | 494 | 0.840 |
| C2HNO_rxn4440_R | NN+RL rl4L | 73.3 | 81.5 | 8.10 (8.10 / 8.10) | 1 | 7.89 ± 0.28 (5) | 97.30 | 608 | 0.853 |
| C2HNO_rxn4440_R | NN+RL rl4n29f | 71.4 | 86.8 | 10.39 (10.39 / 10.39) | 1 | 10.32 ± 0.36 (5) | 96.47 | 477 | 0.878 |
| C2HNO_rxn4440_TS | truncated CCSD | 15.4 | 278.1 | 81.21 (81.21 / 81.21) | 1 | 83.30 ± 2.42 (5) | 74.48 | 80 | 0.970 |
| C2HNO_rxn4440_TS | compressed DF (Lin) | -87.8 | 602.2 | 6.95 (5.84 / 8.02) | 10 | 6.49 ± 0.37 (5) | 98.01 | 930 | 0.426 |
| C2HNO_rxn4440_TS | optimize=True label | 67.1 | 115.8 | 9.99 (9.99 / 9.99) | 1 | 9.40 ± 0.49 (5) | 97.12 | 588 | 0.825 |
| C2HNO_rxn4440_TS | pretrained x4 | 46.2 | 181.2 | 8.83 (8.83 / 8.83) | 1 | 8.24 ± 0.58 (5) | 97.48 | 592 | 0.782 |
| C2HNO_rxn4440_TS | NN+RL rl4L | 68.6 | 110.9 | 7.64 (7.64 / 7.64) | 1 | 7.18 ± 0.50 (5) | 97.80 | 650 | 0.811 |
| C2HNO_rxn4440_TS | NN+RL rl4n29f | 67.9 | 113.0 | 10.83 (10.83 / 10.83) | 1 | 10.79 ± 0.67 (5) | 96.69 | 509 | 0.841 |
| C2HNO_rxn4441_P | truncated CCSD | 13.5 | 238.6 | 85.92 (85.92 / 85.92) | 1 | 89.81 ± 8.64 (5) | 67.26 | 87 | 0.975 |
| C2HNO_rxn4441_P | compressed DF (Lin) | -158.1 | 693.0 | 9.96 (9.39 / 10.50) | 10 | 10.66 ± 0.43 (5) | 96.11 | 1012 | 0.402 |
| C2HNO_rxn4441_P | optimize=True label | 65.8 | 100.2 | 13.12 (13.12 / 13.12) | 1 | 11.72 ± 0.86 (5) | 95.73 | 572 | 0.847 |
| C2HNO_rxn4441_P | pretrained x4 | 51.5 | 138.0 | 12.23 (12.23 / 12.23) | 1 | 12.24 ± 0.37 (5) | 95.54 | 525 | 0.837 |
| C2HNO_rxn4441_P | NN+RL rl4L | 65.7 | 100.3 | 10.17 (10.17 / 10.17) | 1 | 10.27 ± 0.18 (5) | 96.26 | 638 | 0.846 |
| C2HNO_rxn4441_P | NN+RL rl4n29f | 63.5 | 106.0 | 13.02 (13.02 / 13.02) | 1 | 13.10 ± 0.40 (5) | 95.22 | 504 | 0.868 |
| C2HNO_rxn4441_R | truncated CCSD | 17.1 | 243.1 | 75.32 (75.32 / 75.32) | 1 | 81.36 ± 3.47 (5) | 72.18 | 73 | 0.985 |
| C2HNO_rxn4441_R | compressed DF (Lin) | -122.2 | 644.1 | 7.38 (6.75 / 7.99) | 10 | 6.93 ± 0.26 (5) | 97.63 | 1076 | 0.504 |
| C2HNO_rxn4441_R | optimize=True label | 61.6 | 115.2 | 9.16 (9.16 / 9.16) | 1 | 9.53 ± 0.43 (5) | 96.74 | 554 | 0.833 |
| C2HNO_rxn4441_R | pretrained x4 | 56.4 | 130.2 | 8.94 (8.94 / 8.94) | 1 | 9.37 ± 0.30 (5) | 96.80 | 494 | 0.840 |
| C2HNO_rxn4441_R | NN+RL rl4L | 73.3 | 81.5 | 8.10 (8.10 / 8.10) | 1 | 7.89 ± 0.28 (5) | 97.30 | 608 | 0.853 |
| C2HNO_rxn4441_R | NN+RL rl4n29f | 71.4 | 86.8 | 10.39 (10.39 / 10.39) | 1 | 10.32 ± 0.36 (5) | 96.47 | 477 | 0.878 |
| C2HNO_rxn4441_TS | truncated CCSD | 10.0 | 260.5 | 78.38 (78.38 / 78.38) | 1 | 78.31 ± 8.61 (5) | 72.83 | 111 | 0.967 |
| C2HNO_rxn4441_TS | compressed DF (Lin) | -79.6 | 509.5 | 9.85 (9.22 / 11.33) | 10 | 10.22 ± 0.48 (5) | 96.45 | 954 | 0.497 |
| C2HNO_rxn4441_TS | optimize=True label | 46.0 | 160.2 | 12.65 (12.65 / 12.65) | 1 | 12.65 ± 0.53 (5) | 95.61 | 649 | 0.801 |
| C2HNO_rxn4441_TS | pretrained x4 | 46.6 | 158.6 | 12.87 (12.87 / 12.87) | 1 | 12.87 ± 0.28 (5) | 95.53 | 644 | 0.819 |
| C2HNO_rxn4441_TS | NN+RL rl4L | 67.1 | 101.6 | 10.75 (10.75 / 10.75) | 1 | 11.26 ± 0.40 (5) | 96.09 | 741 | 0.837 |
| C2HNO_rxn4441_TS | NN+RL rl4n29f | 67.9 | 99.2 | 15.05 (15.05 / 15.05) | 1 | 15.02 ± 0.35 (5) | 94.79 | 578 | 0.864 |
| C2HNO_rxn4442_P | truncated CCSD | 19.9 | 214.3 | 84.23 (84.23 / 84.23) | 1 | 95.60 ± 7.02 (5) | 64.01 | 83 | 0.973 |
| C2HNO_rxn4442_P | compressed DF (Lin) | -128.8 | 598.6 | 8.49 (8.14 / 9.16) | 10 | 8.62 ± 0.16 (5) | 96.75 | 970 | 0.309 |
| C2HNO_rxn4442_P | optimize=True label | 70.0 | 84.8 | 10.86 (10.86 / 10.86) | 1 | 11.28 ± 0.36 (5) | 95.75 | 596 | 0.863 |
| C2HNO_rxn4442_P | pretrained x4 | 50.0 | 136.4 | 12.11 (12.11 / 12.11) | 1 | 11.56 ± 0.66 (5) | 95.65 | 541 | 0.839 |
| C2HNO_rxn4442_P | NN+RL rl4L | 70.9 | 82.4 | 9.71 (9.71 / 9.71) | 1 | 9.94 ± 0.20 (5) | 96.26 | 634 | 0.861 |
| C2HNO_rxn4442_P | NN+RL rl4n29f | 68.0 | 89.7 | 11.31 (11.31 / 11.31) | 1 | 12.75 ± 1.20 (5) | 95.20 | 522 | 0.878 |
| C2HNO_rxn4442_R | truncated CCSD | 17.1 | 243.1 | 75.32 (75.32 / 75.32) | 1 | 81.36 ± 3.47 (5) | 72.18 | 73 | 0.985 |
| C2HNO_rxn4442_R | compressed DF (Lin) | -122.2 | 644.1 | 7.38 (6.75 / 7.99) | 10 | 6.93 ± 0.26 (5) | 97.63 | 1076 | 0.504 |
| C2HNO_rxn4442_R | optimize=True label | 61.6 | 115.2 | 9.16 (9.16 / 9.16) | 1 | 9.53 ± 0.43 (5) | 96.74 | 554 | 0.833 |
| C2HNO_rxn4442_R | pretrained x4 | 56.4 | 130.2 | 8.94 (8.94 / 8.94) | 1 | 9.37 ± 0.30 (5) | 96.80 | 494 | 0.840 |
| C2HNO_rxn4442_R | NN+RL rl4L | 73.3 | 81.5 | 8.10 (8.10 / 8.10) | 1 | 7.89 ± 0.28 (5) | 97.30 | 608 | 0.853 |
| C2HNO_rxn4442_R | NN+RL rl4n29f | 71.4 | 86.8 | 10.39 (10.39 / 10.39) | 1 | 10.32 ± 0.36 (5) | 96.47 | 477 | 0.878 |
| C2HNO_rxn4442_TS | truncated CCSD | 13.0 | 265.3 | 72.63 (72.63 / 72.63) | 1 | 75.48 ± 2.33 (5) | 75.09 | 132 | 0.947 |
| C2HNO_rxn4442_TS | compressed DF (Lin) | -42.9 | 427.4 | 8.73 (7.91 / 10.03) | 10 | 8.74 ± 0.21 (5) | 97.12 | 1000 | 0.541 |
| C2HNO_rxn4442_TS | optimize=True label | 54.0 | 146.3 | 12.20 (12.20 / 12.20) | 1 | 12.14 ± 0.30 (5) | 95.99 | 733 | 0.793 |
| C2HNO_rxn4442_TS | pretrained x4 | 43.3 | 177.3 | 11.90 (11.90 / 11.90) | 1 | 11.73 ± 0.55 (5) | 96.13 | 711 | 0.792 |
| C2HNO_rxn4442_TS | NN+RL rl4L | 63.9 | 117.8 | 11.98 (11.98 / 11.98) | 1 | 11.32 ± 0.45 (5) | 96.26 | 754 | 0.820 |
| C2HNO_rxn4442_TS | NN+RL rl4n29f | 64.2 | 116.9 | 14.31 (14.31 / 14.31) | 1 | 14.07 ± 1.17 (5) | 95.36 | 622 | 0.844 |

### n17: costs

| stage | median s | n |
|:--|--:|--:|
| t_energy (per state / per draw) | 7.0 | 360 |
| t_sample (per state / per draw) | 8.9 | 360 |
| t_ci (per state / per draw) | 20.3 | 360 |
| diag lin, norb 17 (threads [4]) | 46.9 (max 118) | 300 |
| diag lin, norb 17, cdf_lin (threads [4]) | 137.0 (max 234) | 600 |
| diag n2631g, norb 17 (threads [4]) | 219.2 (max 557) | 60 |
| diag n2631g, norb 17, cdf_lin (threads [4]) | 1063.4 (max 2326) | 12 |
| diag top100, norb 17 (threads [2]) | 7.5 (max 12) | 60 |
| diag top100, norb 17, cdf_lin (threads [2]) | 6.7 (max 9) | 12 |
| diag top300, norb 17 (threads [2]) | 48.6 (max 66) | 60 |
| diag top300, norb 17, cdf_lin (threads [2]) | 45.6 (max 65) | 12 |
