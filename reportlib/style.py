"""Palette and label choices, with the reason, so they are not re-litigated."""

# Estimator families. Grouped by family rather than one colour per arm, because
# the report's claim is about occupancy-vs-Royle-Nichols, not about four
# unrelated methods. Red/blue rather than red/green: red/green is unreadable
# for the ~8% of male readers with deuteranopia.
EST_COLOR = {
    "camera_occ": "#b2182b", "array_occ": "#ef8a62",
    "camera_rn":  "#1b5e9c", "array_rn":  "#67a9cf",
}
EST_LABEL = {
    "camera_occ": "Camera occupancy", "array_occ": "Array occupancy",
    "camera_rn":  "Camera Royle-Nichols", "array_rn": "Array Royle-Nichols",
}
EST_MARKER = {"camera_occ": "o", "array_occ": "s", "camera_rn": "o", "array_rn": "s"}

FOCAL, MUTED = "#1b5e9c", "#b0b7bd"

ABUND_ORDER = ["bobcat_like", "intermediate", "deer_like"]
ABUND_LABEL = {"bobcat_like": "Bobcat-like\n(low)", "intermediate": "Intermediate",
               "deer_like": "Deer-like\n(high)"}

# True generating trend in the estimator sweep. Used to convert the reported
# bias column (estimate - truth) back to an estimate.
TVB_TRUE = -0.1795

SPECIES_ORDER = ["Bobcat", "White-tailed deer", "Moose"]

# Nine covariates: soil_sand is dropped fleet-wide and acts as the implicit
# reference level, because the three soil fractions are compositional and sum
# to a constant, which put a collinearity ridge in moose's posterior.
COVAR_ORDER = ["Human_pop", "NDVI_mean", "Ag", "Deciduous", "Evergreen",
               "Mixed", "terrain_ruggedness", "soil_clay", "soil_silt"]
COVAR_LABEL = {
    "Human_pop": "Human population", "NDVI_mean": "Vegetation greenness (NDVI)",
    "Ag": "Agriculture", "Deciduous": "Deciduous forest",
    "Evergreen": "Evergreen forest", "Mixed": "Mixed forest",
    "terrain_ruggedness": "Terrain ruggedness", "soil_clay": "Soil clay",
    "soil_silt": "Soil silt",
}
