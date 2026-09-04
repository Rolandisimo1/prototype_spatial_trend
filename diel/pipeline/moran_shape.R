## The paper maps diel SHAPE, so residual autocorrelation must be assessed on a shape statistic.
## Array-level mean deviance residual is dominated by abundance (93% zeros; r=-0.82 with total count),
## so it tests the wrong thing. Here we test:
##   (a) the array random harmonic coefficients (the model's own shape deviations) -- if the spatial
##       smooth captured broad-scale shape, these should be spatially unstructured;
##   (b) observed-minus-predicted % nocturnal per array (a direct shape residual).
suppressPackageStartupMessages({library(mgcv); library(data.table); library(spdep)})
source("moran.R")
SUNRISE <- 6.2202; SUNSET <- 17.9998
args <- commandArgs(trailingOnly=TRUE); SP <- args[1]
tag <- gsub("[^A-Za-z0-9]+","_",SP)
MOD <- readRDS("mod_universe.rds"); d <- droplevels(MOD[species==SP]); rm(MOD); gc(FALSE)
d[, `:=`(array_id=factor(array_id), dep_key=factor(dep_key))]
m <- readRDS(sprintf("slim2_%s_ds.rds", tag))

## (a) array random-effect coefficients for the first harmonic
labs <- sapply(m$smooth, function(s) s$label)
cf <- coef(m)
getre <- function(tg) { i <- which(labs==tg); v <- cf[m$smooth[[i]]$first.para:m$smooth[[i]]$last.para]
  data.table(array_id=levels(d$array_id), v=v) }
s1 <- getre("s(sin1,array_id)"); c1 <- getre("s(cos1,array_id)")
xy <- unique(d[, .(array_id=as.character(array_id), X, Y)])[, .(X=mean(X), Y=mean(Y)), by=array_id]
re <- merge(merge(s1, c1, by="array_id", suffixes=c("_sin1","_cos1")), xy, by="array_id")
re[, amp1 := sqrt(v_sin1^2 + v_cos1^2)]
re[, n := 1]

## (b) observed vs predicted % nocturnal per array
night <- function(h) h < SUNRISE | h > SUNSET
re_lab <- function(m) labs[grepl("array_id", labs)]
d[, lp_pop := predict(m, d, type="link", exclude=re_lab(m))]
d[, pred_rate := exp(lp_pop - log_effort_aso)]
obs <- d[, .(obs_night=sum(count[night(bin_center_h)]), obs_tot=sum(count),
             pr_night=sum(pred_rate[night(bin_center_h)]), pr_tot=sum(pred_rate),
             X=mean(X), Y=mean(Y), events=sum(count)), by=array_id]
obs <- obs[obs_tot >= 30]                      # need enough events for a stable observed %
obs[, `:=`(obs_pct=100*obs_night/obs_tot, pred_pct=100*pr_night/pr_tot)]
obs[, resid_pct := obs_pct - pred_pct][, n := 1]

res <- list(sp=SP)
for (nm in c("v_sin1","v_cos1","amp1")) {
  x <- copy(re); setnames(x, nm, "mean_dres")
  mi <- moran_of(x[, .(array_id, mean_dres, X, Y, n)], nsim=999)
  res[[nm]] <- mi
  cat(sprintf("%-22s RE %-7s arrays=%3d I=%+.4f p_mc=%.3g\n", SP, nm, mi$n, mi$I, mi$p_mc))
}
x <- copy(obs); setnames(x, "resid_pct", "mean_dres")
mi <- moran_of(x[, .(array_id, mean_dres, X, Y, n)], nsim=999)
res$resid_pct <- mi
cat(sprintf("%-22s obs-pred %%noct arrays=%3d I=%+.4f p_mc=%.3g | sd(resid)=%.2f pt  cor(obs,pred)=%.3f\n",
    SP, mi$n, mi$I, mi$p_mc, sd(obs$resid_pct), cor(obs$obs_pct, obs$pred_pct)))
res$obs <- obs; res$re <- re
saveRDS(res, sprintf("shape_moran_%s.rds", tag))
