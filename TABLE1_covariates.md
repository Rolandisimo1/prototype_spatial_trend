# Table 1. Covariates

| Covariate | Role | Spatial scale | Temporal | Transform | Source dataset | Reference |
|---|---|---|---|---|---|---|
| `Human_pop` | intensity + occupancy | cell50 (50 km) | annual | <bound method Series.transform of covariate                          Human_pop
role                   intensity + occupancy
spatial_scale                 cell50 (50 km)
temporal_resolution                   annual
transform                       standardised
source_dataset                              
reference                                   
status                         SOURCE NEEDED
Name: 0, dtype: str> | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `NDVI_mean` | intensity + occupancy | cell50 (50 km) | annual | <bound method Series.transform of covariate                          NDVI_mean
role                   intensity + occupancy
spatial_scale                 cell50 (50 km)
temporal_resolution                   annual
transform                       standardised
source_dataset                              
reference                                   
status                         SOURCE NEEDED
Name: 1, dtype: str> | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `Ag` | intensity + occupancy | cell50 (50 km) | annual | <bound method Series.transform of covariate                                 Ag
role                   intensity + occupancy
spatial_scale                 cell50 (50 km)
temporal_resolution                   annual
transform                       standardised
source_dataset                              
reference                                   
status                         SOURCE NEEDED
Name: 2, dtype: str> | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `Deciduous` | intensity + occupancy | cell50 (50 km) | annual | <bound method Series.transform of covariate                          Deciduous
role                   intensity + occupancy
spatial_scale                 cell50 (50 km)
temporal_resolution                   annual
transform                       standardised
source_dataset                              
reference                                   
status                         SOURCE NEEDED
Name: 3, dtype: str> | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `Evergreen` | intensity + occupancy | cell50 (50 km) | annual | <bound method Series.transform of covariate                          Evergreen
role                   intensity + occupancy
spatial_scale                 cell50 (50 km)
temporal_resolution                   annual
transform                       standardised
source_dataset                              
reference                                   
status                         SOURCE NEEDED
Name: 4, dtype: str> | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `Mixed` | intensity + occupancy | cell50 (50 km) | annual | <bound method Series.transform of covariate                              Mixed
role                   intensity + occupancy
spatial_scale                 cell50 (50 km)
temporal_resolution                   annual
transform                       standardised
source_dataset                              
reference                                   
status                         SOURCE NEEDED
Name: 5, dtype: str> | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `terrain_ruggedness` | intensity + occupancy | cell50 (50 km) | static | <bound method Series.transform of covariate                 terrain_ruggedness
role                   intensity + occupancy
spatial_scale                 cell50 (50 km)
temporal_resolution                   static
transform                       standardised
source_dataset                              
reference                                   
status                         SOURCE NEEDED
Name: 6, dtype: str> | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `soil_clay` | intensity + occupancy | cell50 (50 km) | static | <bound method Series.transform of covariate                          soil_clay
role                   intensity + occupancy
spatial_scale                 cell50 (50 km)
temporal_resolution                   static
transform                       standardised
source_dataset                              
reference                                   
status                         SOURCE NEEDED
Name: 7, dtype: str> | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `soil_silt` | intensity + occupancy | cell50 (50 km) | static | <bound method Series.transform of covariate                          soil_silt
role                   intensity + occupancy
spatial_scale                 cell50 (50 km)
temporal_resolution                   static
transform                       standardised
source_dataset                              
reference                                   
status                         SOURCE NEEDED
Name: 8, dtype: str> | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `MWMT` | spatially varying coefficient | cell100 field, cell50 data | annual | <bound method Series.transform of covariate                                       MWMT
role                   spatially varying coefficient
spatial_scale             cell100 field, cell50 data
temporal_resolution                           annual
transform                               standardised
source_dataset                                      
reference                                           
status                                 SOURCE NEEDED
Name: 9, dtype: str> | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `MCMT` | spatially varying coefficient | cell100 field, cell50 data | annual | <bound method Series.transform of covariate                                       MCMT
role                   spatially varying coefficient
spatial_scale             cell100 field, cell50 data
temporal_resolution                           annual
transform                               standardised
source_dataset                                      
reference                                           
status                                 SOURCE NEEDED
Name: 10, dtype: str> | **[SOURCE NEEDED]** | **[REFERENCE NEEDED]** |
| `yday` | detection | per deployment-occasion | occasion | <bound method Series.transform of covariate                                            yday
role                                            detection
spatial_scale                     per deployment-occasion
temporal_resolution                              occasion
transform              standardised; linear and quadratic
source_dataset                 camera deployment metadata
reference                   derived from deployment dates
status                                         DOCUMENTED
Name: 11, dtype: str> | camera deployment metadata | derived from deployment dates |
| `canopy_height` | detection | per deployment | static | <bound method Series.transform of covariate                                                                  canopy_height
role                                                                           detection
spatial_scale                                                             per deployment
temporal_resolution                                                               static
transform                                                                   standardised
source_dataset                                                      canopy height raster
reference              cited for the same camera dataset by Goldstein et al. (in review)
status                                                              CONFIRM WITH ARIELLE
Name: 12, dtype: str> | canopy height raster | cited for the same camera dataset by Goldstein et al. (in review) |
| `log_roaddist` | detection | per deployment | static | <bound method Series.transform of covariate                                               log_roaddist
role                                                       detection
spatial_scale                                         per deployment
temporal_resolution                                           static
transform                                      log then standardised
source_dataset                 Global Roads Inventory Project (GRIP)
reference              named in rebuild_detections.R; version needed
status                                                VERSION NEEDED
Name: 13, dtype: str> | Global Roads Inventory Project (GRIP) | named in rebuild_detections.R; version needed |

Covariates removed in the two documented reductions from the initial set of seventeen: `CWD`, `NDVI_sd`, `Impervious`, `PDSI`, `PPT`, `Temp`, `elevation`, `soil_sand`.
