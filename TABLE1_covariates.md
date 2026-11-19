# Table 1. Covariates

| Covariate | Role | Spatial scale | Temporal | Transform | Source dataset | Reference |
|---|---|---|---|---|---|---|
| `Human_pop` | intensity + occupancy | cell50 (50 km) | annual | standardised | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `NDVI_mean` | intensity + occupancy | cell50 (50 km) | annual | standardised | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `Ag` | intensity + occupancy | cell50 (50 km) | annual | standardised | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `Deciduous` | intensity + occupancy | cell50 (50 km) | annual | standardised | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `Evergreen` | intensity + occupancy | cell50 (50 km) | annual | standardised | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `Mixed` | intensity + occupancy | cell50 (50 km) | annual | standardised | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `terrain_ruggedness` | intensity + occupancy | cell50 (50 km) | static | standardised | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `soil_clay` | intensity + occupancy | cell50 (50 km) | static | standardised | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `soil_silt` | intensity + occupancy | cell50 (50 km) | static | standardised | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `MWMT` | spatially varying coefficient | cell100 field, cell50 data | annual | standardised | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `MCMT` | spatially varying coefficient | cell100 field, cell50 data | annual | standardised | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `yday` | detection | per deployment-occasion | occasion | standardised; linear and quadratic | camera deployment metadata | derived from deployment dates |
| `canopy_height` | detection | per deployment | static | standardised | canopy height raster | cited for the same camera dataset by Goldstein et al. (in review) |
| `log_roaddist` | detection | per deployment | static | log then standardised | Global Roads Inventory Project (GRIP) | named in rebuild_detections.R; version needed |
