## theta0 / theta1 — all 8 models

| model | parameter | mean | median | sd | q025 | q975 | R-hat |
|---|---|---|---|---|---|---|---|
| moose_v1fix9 eco | `theta0` | -3.8856 | -3.8795 | 0.1483 | -4.1987 | -3.6094 | 1.020 |
| moose_v1fix9 ns | `theta0` | -3.9142 | -3.9061 | 0.1570 | -4.2443 | -3.6276 | 1.011 |
| moose_v2b eco | `theta0` | -3.9539 | -3.9507 | 0.1633 | -4.2891 | -3.6417 | 1.030 |
| moose_v2b ns | `theta0` | -4.0041 | -3.9945 | 0.1834 | -4.3823 | -3.6685 | 1.030 |
| WTD_v2b eco | `theta0` | -4.0997 | -4.1030 | 0.0842 | -4.2586 | -3.9240 | 1.018 |
| WTD_v2b ns | `theta0` | -4.1079 | -4.1120 | 0.0853 | -4.2727 | -3.9370 | 1.039 |
| bobcat_v2b ns | `theta0` | -6.0183 | -6.0168 | 0.1021 | -6.2196 | -5.8214 | 1.006 |
| bobcat_v2b eco | `theta0` | -6.0294 | -6.0279 | 0.1030 | -6.2345 | -5.8346 | 1.006 |
| bobcat_v2b eco | `theta1` | 0.6146 | 0.6140 | 0.0268 | 0.5640 | 0.6689 | 1.007 |
| bobcat_v2b ns | `theta1` | 0.6127 | 0.6123 | 0.0268 | 0.5605 | 0.6658 | 1.006 |
| moose_v2b ns | `theta1` | 0.5285 | 0.5249 | 0.0498 | 0.4405 | 0.6313 | 1.024 |
| moose_v2b eco | `theta1` | 0.5256 | 0.5253 | 0.0435 | 0.4432 | 0.6120 | 1.036 |
| moose_v1fix9 ns | `theta1` | 0.4713 | 0.4682 | 0.0457 | 0.3882 | 0.5676 | 1.014 |
| moose_v1fix9 eco | `theta1` | 0.4662 | 0.4660 | 0.0428 | 0.3859 | 0.5565 | 1.068 |
| WTD_v2b ns | `theta1` | 0.4115 | 0.4123 | 0.0158 | 0.3797 | 0.4414 | 1.041 |
| WTD_v2b eco | `theta1` | 0.4098 | 0.4103 | 0.0157 | 0.3772 | 0.4396 | 1.019 |

## occ_beta — 9 covariates x 8 models

### bobcat_v2b ns

| # | covariate | mean | median | sd | q025 | q975 | R-hat |
|---|---|---|---|---|---|---|---|
| 1 | Human_pop | -0.3796 | -0.3789 | 0.0437 | -0.4671 | -0.2954 | 1.001 |
| 2 | NDVI_mean | 0.1171 | 0.1172 | 0.0564 | 0.0054 | 0.2257 | 1.003 |
| 3 | Ag | 0.0483 | 0.0483 | 0.0303 | -0.0107 | 0.1084 | 1.002 |
| 4 | Deciduous | 0.2170 | 0.2166 | 0.0599 | 0.1010 | 0.3361 | 1.001 |
| 5 | Evergreen | -0.0065 | -0.0071 | 0.0476 | -0.0987 | 0.0880 | 1.002 |
| 6 | Mixed | 0.0441 | 0.0438 | 0.0387 | -0.0307 | 0.1204 | 1.000 |
| 7 | terrain_ruggedness | 0.0739 | 0.0738 | 0.0223 | 0.0307 | 0.1181 | 1.001 |
| 8 | soil_clay | -0.1265 | -0.1263 | 0.0385 | -0.2020 | -0.0517 | 1.001 |
| 9 | soil_silt | -0.0130 | -0.0127 | 0.0416 | -0.0950 | 0.0691 | 1.002 |

### bobcat_v2b eco

| # | covariate | mean | median | sd | q025 | q975 | R-hat |
|---|---|---|---|---|---|---|---|
| 1 | Human_pop | -0.3798 | -0.3795 | 0.0432 | -0.4652 | -0.2961 | 1.001 |
| 2 | NDVI_mean | 0.1105 | 0.1097 | 0.0548 | 0.0044 | 0.2198 | 1.003 |
| 3 | Ag | 0.0591 | 0.0590 | 0.0303 | -0.0004 | 0.1190 | 1.001 |
| 4 | Deciduous | 0.2195 | 0.2202 | 0.0614 | 0.0985 | 0.3409 | 1.002 |
| 5 | Evergreen | -0.0112 | -0.0111 | 0.0477 | -0.1041 | 0.0828 | 1.002 |
| 6 | Mixed | 0.0445 | 0.0446 | 0.0387 | -0.0311 | 0.1206 | 1.001 |
| 7 | terrain_ruggedness | 0.0690 | 0.0689 | 0.0222 | 0.0259 | 0.1127 | 1.000 |
| 8 | soil_clay | -0.0898 | -0.0895 | 0.0395 | -0.1681 | -0.0124 | 1.001 |
| 9 | soil_silt | -0.0337 | -0.0339 | 0.0422 | -0.1163 | 0.0486 | 1.002 |

### WTD_v2b ns

| # | covariate | mean | median | sd | q025 | q975 | R-hat |
|---|---|---|---|---|---|---|---|
| 1 | Human_pop | -0.5013 | -0.5010 | 0.0297 | -0.5603 | -0.4442 | 1.000 |
| 2 | NDVI_mean | -0.0448 | -0.0451 | 0.0340 | -0.1110 | 0.0232 | 1.003 |
| 3 | Ag | 0.0111 | 0.0110 | 0.0208 | -0.0293 | 0.0523 | 1.001 |
| 4 | Deciduous | 0.0655 | 0.0653 | 0.0351 | -0.0031 | 0.1345 | 1.000 |
| 5 | Evergreen | 0.1431 | 0.1435 | 0.0310 | 0.0824 | 0.2037 | 1.002 |
| 6 | Mixed | 0.0020 | 0.0020 | 0.0214 | -0.0403 | 0.0439 | 1.002 |
| 7 | terrain_ruggedness | -0.0433 | -0.0433 | 0.0162 | -0.0750 | -0.0112 | 1.001 |
| 8 | soil_clay | 0.0352 | 0.0351 | 0.0245 | -0.0127 | 0.0836 | 1.001 |
| 9 | soil_silt | 0.0330 | 0.0328 | 0.0286 | -0.0220 | 0.0897 | 1.001 |

### WTD_v2b eco

| # | covariate | mean | median | sd | q025 | q975 | R-hat |
|---|---|---|---|---|---|---|---|
| 1 | Human_pop | -0.4994 | -0.4989 | 0.0302 | -0.5599 | -0.4408 | 1.001 |
| 2 | NDVI_mean | -0.0445 | -0.0446 | 0.0335 | -0.1105 | 0.0220 | 1.003 |
| 3 | Ag | 0.0130 | 0.0129 | 0.0211 | -0.0283 | 0.0539 | 1.002 |
| 4 | Deciduous | 0.0729 | 0.0728 | 0.0358 | 0.0029 | 0.1442 | 1.001 |
| 5 | Evergreen | 0.1448 | 0.1448 | 0.0309 | 0.0848 | 0.2060 | 1.000 |
| 6 | Mixed | 0.0046 | 0.0044 | 0.0215 | -0.0372 | 0.0471 | 1.001 |
| 7 | terrain_ruggedness | -0.0442 | -0.0441 | 0.0164 | -0.0764 | -0.0121 | 1.000 |
| 8 | soil_clay | 0.0438 | 0.0436 | 0.0253 | -0.0053 | 0.0941 | 1.000 |
| 9 | soil_silt | 0.0298 | 0.0299 | 0.0295 | -0.0284 | 0.0879 | 1.004 |

### moose_v2b ns

| # | covariate | mean | median | sd | q025 | q975 | R-hat |
|---|---|---|---|---|---|---|---|
| 1 | Human_pop | -3.7734 | -3.7078 | 0.9479 | -5.7953 | -2.0903 | 1.016 |
| 2 | NDVI_mean | -0.4943 | -0.4943 | 0.1073 | -0.7054 | -0.2866 | 1.021 |
| 3 | Ag | 0.1131 | 0.1130 | 0.0614 | -0.0076 | 0.2341 | 1.001 |
| 4 | Deciduous | 0.0626 | 0.0639 | 0.1958 | -0.3218 | 0.4491 | 1.009 |
| 5 | Evergreen | 1.0210 | 1.0187 | 0.1943 | 0.6402 | 1.4133 | 1.007 |
| 6 | Mixed | 0.5104 | 0.5069 | 0.1745 | 0.1748 | 0.8662 | 1.002 |
| 7 | terrain_ruggedness | -0.1382 | -0.1389 | 0.0737 | -0.2820 | 0.0069 | 1.000 |
| 8 | soil_clay | -0.0487 | -0.0483 | 0.1625 | -0.3668 | 0.2699 | 1.004 |
| 9 | soil_silt | -0.3929 | -0.3818 | 0.2727 | -0.9583 | 0.1115 | 1.001 |

### moose_v2b eco

| # | covariate | mean | median | sd | q025 | q975 | R-hat |
|---|---|---|---|---|---|---|---|
| 1 | Human_pop | -3.7978 | -3.7588 | 0.9146 | -5.7210 | -2.0955 | 1.003 |
| 2 | NDVI_mean | -0.5108 | -0.5109 | 0.1038 | -0.7186 | -0.3070 | 1.004 |
| 3 | Ag | 0.1124 | 0.1129 | 0.0601 | -0.0072 | 0.2292 | 1.006 |
| 4 | Deciduous | 0.0736 | 0.0769 | 0.1905 | -0.3077 | 0.4324 | 1.003 |
| 5 | Evergreen | 1.0588 | 1.0553 | 0.1936 | 0.6970 | 1.4556 | 1.006 |
| 6 | Mixed | 0.5425 | 0.5429 | 0.1813 | 0.1847 | 0.8969 | 1.024 |
| 7 | terrain_ruggedness | -0.1408 | -0.1409 | 0.0736 | -0.2848 | 0.0039 | 1.000 |
| 8 | soil_clay | -0.0972 | -0.0973 | 0.1627 | -0.4190 | 0.2230 | 1.003 |
| 9 | soil_silt | -0.4198 | -0.4118 | 0.2794 | -0.9974 | 0.0969 | 1.004 |

### moose_v1fix9 ns

| # | covariate | mean | median | sd | q025 | q975 | R-hat |
|---|---|---|---|---|---|---|---|
| 1 | Human_pop | -2.1125 | -2.0605 | 0.8049 | -3.8239 | -0.6852 | 1.001 |
| 2 | NDVI_mean | -0.2209 | -0.2247 | 0.1364 | -0.4773 | 0.0585 | 1.004 |
| 3 | Ag | 0.1155 | 0.1159 | 0.0673 | -0.0192 | 0.2461 | 1.005 |
| 4 | Deciduous | -0.3103 | -0.3000 | 0.2655 | -0.8708 | 0.1831 | 1.005 |
| 5 | Evergreen | 0.6103 | 0.6130 | 0.2476 | 0.1053 | 1.0820 | 1.011 |
| 6 | Mixed | 0.0931 | 0.0935 | 0.2284 | -0.3543 | 0.5427 | 1.001 |
| 7 | terrain_ruggedness | -0.1240 | -0.1244 | 0.0857 | -0.2914 | 0.0442 | 1.001 |
| 8 | soil_clay | -0.0117 | -0.0119 | 0.1844 | -0.3730 | 0.3486 | 1.006 |
| 9 | soil_silt | -0.2906 | -0.2797 | 0.2827 | -0.8742 | 0.2385 | 1.000 |

### moose_v1fix9 eco

| # | covariate | mean | median | sd | q025 | q975 | R-hat |
|---|---|---|---|---|---|---|---|
| 1 | Human_pop | -2.2580 | -2.2082 | 0.8477 | -4.0703 | -0.7381 | 1.005 |
| 2 | NDVI_mean | -0.2270 | -0.2285 | 0.1385 | -0.4981 | 0.0478 | 1.002 |
| 3 | Ag | 0.1151 | 0.1160 | 0.0695 | -0.0232 | 0.2510 | 1.001 |
| 4 | Deciduous | -0.3459 | -0.3389 | 0.2735 | -0.8985 | 0.1628 | 1.018 |
| 5 | Evergreen | 0.6307 | 0.6300 | 0.2499 | 0.1469 | 1.1243 | 1.003 |
| 6 | Mixed | 0.1030 | 0.1041 | 0.2336 | -0.3529 | 0.5608 | 1.003 |
| 7 | terrain_ruggedness | -0.1241 | -0.1243 | 0.0854 | -0.2921 | 0.0440 | 1.001 |
| 8 | soil_clay | -0.0365 | -0.0386 | 0.1868 | -0.3990 | 0.3318 | 1.003 |
| 9 | soil_silt | -0.3007 | -0.2935 | 0.2875 | -0.8834 | 0.2388 | 1.002 |


---

## theta1 vs TRI — the ranking disagreement

| species (ns fit) | theta1 | TRI | rank theta1 | rank TRI |
|---|---|---|---|---|
| bobcat_v2b | 0.6127 | 0.0005 | 1 | 4 |
| moose_v2b | 0.5285 | 0.2550 | 2 | 3 |
| moose_v1fix9 | 0.4713 | 0.4880 | 3 | 2 |
| WTD_v2b | 0.4115 | 0.7836 | 4 | 1 |

**Spearman = -1.000, Pearson = -0.993.** The ranks are exactly reversed.

The working hypothesis -- TRI measures temporal agreement of the two streams
over 18 years, theta1 measures cross-sectional scaling of iNat intensity
against camera-estimated abundance -- is sound as a statement about what the
two quantities ARE, and it correctly implies a species can score high on one
and low on the other without contradiction.

But it predicts only that they CAN diverge. It does not predict a perfect
inversion, and that is what the data show. Distinctness alone would produce
a scatter, not a monotone reversal across all four fits.

One natural mechanism was tested and REFUTED: that theta1 is attenuated
toward zero in species with little spatial contrast in abundance (a
regression-dilution argument -- flat predictor, shallow slope). Spatial
dispersion of filtered records across occupied cell50:

| species | occupied cells | records | sd(log count) | IQR(log) | Gini | max/median |
|---|---|---|---|---|---|---|
| bobcat | 1,390 | 18,725 | 1.322 | 1.946 | 0.790 | 414 |
| moose | 421 | 11,581 | 1.531 | 2.197 | 0.785 | 259 |
| WTD | 2,104 | 188,939 | **1.789** | 2.593 | 0.797 | 633 |

The relationship runs OPPOSITE to the attenuation prediction: WTD has the
LARGEST spatial contrast and the LOWEST theta1; bobcat the smallest contrast
and the highest theta1. So low dynamic range does not explain it. theta1 is
also not monotone in record count or in occupied-cell count (moose has the
fewest of both yet sits mid-range on theta1).

RECOMMENDATION for the report: state the conceptual distinction, which is
correct and worth making, and report the inversion as an observation that is
NOT explained -- rather than presenting the distinction as if it accounted
for it. With three species / four fits, Spearman = -1 carries p ~ 0.083
two-sided; it is entirely possible this is coincidence. Asserting a
mechanism here would be over-reading n=4.
