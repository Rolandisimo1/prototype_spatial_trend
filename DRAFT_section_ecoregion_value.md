# Is the ecoregion model worth fitting?

We fit every species twice, under two spatial structures for the population
trend. Both use the same data, the same nine habitat covariates, and the same
integration of camera detections with iNaturalist records. They differ in one
respect: how the trend is allowed to vary across the country.

The **national model** estimates a single trend and applies it everywhere. The
**ecoregion model** adds one deviation per EPA Level I ecoregion, drawn from a
shared normal distribution whose width (`sigma_region`) is itself estimated.
Regions with little data are pulled toward the national mean; regions with
plenty of data are allowed to depart from it. This is standard partial pooling,
and the national trend parameter is defined identically in both models
(`total_var_beta = year_beta + year_var`), with the regional deviations centred
on zero. The two models therefore estimate the same national quantity, which
makes them directly comparable.

That comparison is worth making explicitly, because the regional layer is not
free.

## The cost is a real loss of national precision

Adding eight regional deviations widens the 95% interval on the national trend
by a factor of 1.9 to 3.1:

| Species | National model | Ecoregion model | Interval width |
|---|---|---|---|
| Bobcat | −0.242 (−0.314, −0.172) | −0.145 (−0.372, +0.066) | 3.1× wider |
| White-tailed deer | +0.177 (+0.136, +0.221) | +0.187 (+0.109, +0.266) | 1.9× wider |
| Moose | −0.121 (−0.227, −0.025) | −0.118 (−0.342, +0.088) | 2.1× wider |

For bobcat and moose this changes the conclusion, not just its precision. Under
the national model both species show a decline whose interval excludes zero.
Under the ecoregion model neither does. A reader given only the ecoregion fit
would conclude that we cannot establish a national trend for either species.

The two fits do not contradict each other. Every ecoregion interval contains the
corresponding national estimate, and all three pairs of intervals overlap. The
national mean trend simply becomes harder to pin down once the model must also
decide how much of the signal belongs to regional variation, because the mean
and the deviations are estimated jointly and trade off against one another.

## Whether the cost buys anything depends on the species

The regional layer is worth its cost only if the regional estimates say
something the national estimate cannot. Three measurements settle that, and they
give three different answers.

**Is there regional structure to resolve?** `sigma_region` is the direct readout.
For bobcat it is 0.294 with a 95% interval of (0.150, 0.557) — comfortably away
from zero, so regions genuinely differ. For white-tailed deer it is 0.111
(0.010, 0.322) and for moose 0.202 (0.026, 0.533); both intervals run down
almost to zero, meaning the data are consistent with little or no regional
variation.

**Are the regional estimates distinguishable from one another?** Counting only
regions whose own interval excludes zero, and then counting pairs of those
regions whose intervals do not overlap: bobcat has 6 clear regions and 5 of 15
pairs distinguishable. Deer has 5 clear regions and **zero** distinguishable
pairs. Moose has 1 clear region, so no pairs are possible. For deer, in other
words, the model produces eight regional numbers that are all saying the same
thing as each other and as the national estimate.

**Does any region contradict the national picture?** This is the finding a
national model structurally cannot produce, and it occurs once. For bobcat,
Northern Forests shows an increase of +0.310 (+0.105, +0.519) against a national
decline of −0.242. That region carries 1,450 camera sites and 101 iNaturalist
cells, so the estimate is informed by real local data rather than produced by
shrinkage.

## Recommendation

**Report the national model as the headline for all three species, and the
ecoregion model only where it earns its place.**

- **Bobcat: fit and report both.** The regional layer pays for itself. Bobcat is
  not declining uniformly — it is declining across most of its range while
  increasing in the Northern Forests, and that reversal is invisible to the
  national model. The wider national interval is the price of learning it, and
  it is worth paying.
- **White-tailed deer: report the national model.** The regional layer confirms
  spatial uniformity, which is a genuine result and worth one sentence, but it
  yields no region distinguishable from any other. Deer is the species where the
  precision is spent for nothing.
- **Moose: report the national model.** At moose's data volume (1,654 camera
  sites against roughly 20,000 for the other two species) the regional layer
  cannot resolve regional structure. Only one of eight regions supports a
  direction, and the national conclusion is lost in the process.

## Two caveats on reading regional estimates

**Regions with no data are not evidence.** Partial pooling produces an estimate
for every region whether or not that region contains data. Mediterranean
California has no camera sites and no iNaturalist cells for deer or moose, and
Marine West Coast Forest has none for moose. Those estimates are the national
mean plus shrinkage, and should not be read as statements about those regions.
They are greyed out in the regional trend maps.

**Fitting both models is a defensible diagnostic, not duplicated effort.** You
cannot know in advance which of the three cases above you are in, and
`sigma_region` is what tells you. The recommendation here is about what to
*report*, not about what to fit. Running both and reading `sigma_region` before
choosing is the sound procedure, and is what we did.
