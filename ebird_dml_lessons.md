# Lessons from the eBird Double-ML trend model for our integrated SVC trend work

**Source:** Fink et al., *A Double Machine Learning Trend Model for Citizen Science Data*, Methods Ecol. Evol. 2023 (doi:10.1111/2041-210X.14186); code `trendswf` v1.0.0, Zenodo 8092408.
**For:** Roland / Krishna / Arielle — deciding whether anything from the eBird approach should shape the ecoregion (EPA Level I) spatially-varying-trend prototype.

---

## 1. What the eBird DML model is (so we know what's borrowable)

- **Two-stage hurdle**: an occurrence-rate model + a count-given-occurrence model, combined via the chain rule into an abundance percent-per-year change. (Loosely parallels our camera-occupancy + iNat-intensity split.)
- **The trend is a causal-forest treatment effect.** Year = "treatment." The trend `τ(W)` is a *nonparametric function of covariates W*, estimated by a random forest. Spatial variation in the trend emerges because W contains spatial/habitat features and the forest lets the year-effect vary across them — **no grid, no CAR field, no regions.** Spatial resolution is entirely covariate-driven.
- **Interannual effort change is handled as a confounder**, three ways:
  1. a **propensity model** (predict survey year from habitat + effort features),
  2. a fixed **"prediction unit"** — predict with effort covariates held at standardized reference values, so mapped trend is at constant effort,
  3. a **residual-confounding diagnostic** (their headline contribution): simulate a known-trend null with matched confounding, measure leftover bias, subtract it.

## 2. Top-line recommendation: borrow lessons, not the method

Causal forests need eBird's enormous checklist volume to fit a nonparametric trend surface without overfitting. Our camera + iNat mammal data is far sparser. Switching to DML would discard the two things that make our model valuable — **Bayesian fusion of cameras + iNaturalist**, and **full uncertainty propagation** through the CAR structure. **Do not re-architect toward DML.** Extract the four lessons below.

## 3. Four transferable lessons (all cheap, none a redesign)

1. **Add a spatially-CONSTANT (null) scenario to the simulation.** The paper validates that the method distinguishes *spatially constant from spatially varying* trends. Our current sim tests "can it recover real structure"; add the complementary test: simulate a truth where every ecoregion trends identically, fit the ecoregion model, confirm it does **not** invent regional structure (region deviations → ~0, `sigma_region` small). Guards against false-positive spatial structure — reviewers will ask.
2. **Map trends at a fixed reference effort — verify we already do.** Their "prediction unit" fixes effort at prediction time so the map shows population signal, not effort artifact. Our `mu` (latent intensity, with `inat_effort` in the `dnbinom` denominator) is *already* the effort-standardized quantity, so we're structurally fine — but confirm the prediction code maps `mu` at a fixed effort, not at each year's observed (rising) effort.
3. **Residual-confounding diagnostic — optional validation, not a model change.** Simulate iNat data with a known interannual effort trend but ZERO population trend; confirm the 50km offset absorbs it and the estimated trend stays ~0. A "prove the offset works" receipt, if/when wanted. (Not reopening the offset decision — this only documents that it holds.)
4. **Foreground the advantage eBird does NOT have.** Their entire methodological burden exists *because citizen science is all they have*. We fuse **camera traps — an effort-controlled gold standard — into the same model.** That is a stronger design for precisely the problem their paper addresses, and it is the lead argument for why our integrated approach is defensible against the "citizen-science effort is changing" critique.

## 4. On mechanism of spatial flexibility (context, not an action)

eBird's trend heterogeneity is covariate-driven (`τ(W)`), not geographically partitioned. Our EPA Level I ecoregion random effect is a coarse, structured version of the same idea. A year×covariate interaction would give finer, higher-power resolution (one global coefficient, no per-cell shrinkage) — noted here for completeness; deliberately deferred to keep this prototype's simulation simple.

## 5. Net effect on the prototype plan

Unchanged core: EPA Level I ecoregion random-effect deviation on the trend, validated by simulation before any refit. **One concrete addition** (Lesson 1: null/constant-trend scenario) and **one confirmation** (Lesson 2: predict at fixed effort). Lessons 3–4 are framing/validation for the Krishna+Arielle conversation, not code.
