# One-name-per-species filter in `01_Prep_UMFs.R` drops renamed camera records

**Draft for Roland to send. Not sent.** 2026-10-08

Brittany, this is a heads-up and the evidence behind it. We run on a corrected copy in the meantime, so nothing is blocked on our side.

## What happens

`Data_prep/01_Prep_UMFs.R`, line 150:

```r
sequences <- sequences[sciname %in% species_list]
```

`species_list` comes from `species_list.csv`, which has one scientific name per species. When Wildlife Insights changes a species' name partway through the record, every sequence filed under the other name is dropped before the detection histories are built. No error or warning is raised. The same filter dropped sequences carrying a wrong name.

## Evidence from `combined_sequences_all.csv`

| Species | Name in `species_list.csv` | Other name in the data | Sequences (deployments) under the other name | Years |
|---|---|---|---|---|
| Collared peccary | *Pecari tajacu*: 1,640 (151), 2015–2023 | *Dicotyles tajacu* | 1,510 (123) | 2024–2025 |
| Long-tailed weasel | *Mustela frenata*: 184 (105), 2012–2023 | *Neogale frenata* | 186 (93) | 2024–2025 |
| Fisher | *Pekania pennanti*: 545 (319) | *Martes pennanti* | 2 (2) | 2013 |
| Elk | *Cervus canadensis*: 11,540 (883) | *Cervus elaphus*, a labelling error | 855 (66) | 2017–2021 |

- **Peccary:** the name changed cleanly at 2023/24, with no overlap between the two names. Keyed on *Pecari tajacu*, the record loses **47.9%** of its sequences, all of them at the recent end. In `umflist.RDS` peccary has 0 detections in 2024 (0 of 3,059 deployments) and in 2025 (0 of 3,520), after 62 detecting deployments in 2023. In our trend model that showed up as a **−54%/yr decline**. With both names, 2024 has 41 detecting deployments and 2025 has 74.
- **Long-tailed weasel:** the same clean 2023/24 break, *Mustela* → *Neogale*. Half the record (50.3%) sits under the new name.
- **Elk:** *C. elaphus* and *C. canadensis* are separate species in the Mammal Diversity Database (MDD v2.5), so no name lookup will merge them. No wild *C. elaphus* occurs in North America, and all 855 sequences are inside North America. We reassign them by a recorded decision; the *Cervus nippon* (sika deer) on Long Island is real and stays separate.
- **Feral hog:** besides *Sus scrofa*, 20 sequences are filed as *Sus scrofa scrofa*.
- More renames will come as WI follows MDD. Across the whole file, 485 distinct camera names don't resolve to any MDD species.

## How we checked our correction

We rebuilt the detection matrices downstream with your window logic unchanged (10–100 inner days, full 10-day windows), keeping your deployments, rows and covariates. As a gate, rebuilding every species from its single original name reproduces `umflist.RDS` exactly for all 36 species: zero differing cells and identical NA patterns. Only then were the extra names added, for peccary and elk. That gives us confidence the only difference is the names.

A cheap check that catches this whatever the cause: compare each species' detecting deployments in the camera record with its detecting deployments in the raw sequences, inside CONUS. The fleet median is about 0.90, because of the window filters. Peccary was 0.52, the only outlier; with both names it's 0.94.

## Suggested fix

Key `species_list` on an MDD species id with a **list** of names per species, ideally with the years each name was used. Then union over that list at line 150, i.e. resolve each name to its id and group by id instead of comparing name strings. Two names that resolve to two different MDD species, like elk, need a recorded decision rather than a lookup.

We're glad to share the rebuild script, the name map and the check if they help.
