## Phase 1 track 1 -- FINAL fits: ti() decomposition, nb(), array harmonic REs, sun-time offset.
## Run one species+basis per process so each ~450 MB bam object is freed at exit.
## Usage: Rscript fit_final.R "<species>" <basis ds|gp> [gp_range_km] [extra_term]
suppressPackageStartupMessages({library(mgcv); library(data.table)})
SUNRISE <- 6.2202; SUNSET <- 17.9998; H <- seq(0.25, 23.75, by=0.5)
args <- commandArgs(trailingOnly=TRUE)
SP <- args[1]; BASIS <- args[2]
GPR <- if (length(args) >= 3 && nzchar(args[3])) as.numeric(args[3]) else NA
EXTRA <- if (length(args) >= 4 && nzchar(args[4])) args[4] else NULL
TAG <- paste0(gsub("[^A-Za-z0-9]+","_",SP), "_", BASIS,
              if (!is.na(GPR)) paste0("_", GPR) else "",
              if (!is.null(EXTRA)) "_eco" else "")

MOD <- readRDS("mod_universe.rds")
d <- droplevels(MOD[species==SP]); rm(MOD); gc(FALSE)
if (!is.null(EXTRA)) {
  eco <- fread("dep_ecoregion_l1.csv")
  d <- merge(d, eco, by="dep_key", all.x=TRUE)
  d[, ecoregion_l1 := factor(ecoregion_l1)]
  stopifnot(!any(is.na(d$ecoregion_l1)))
}
d[, `:=`(array_id=factor(array_id), dep_key=factor(dep_key))]

RE <- paste("s(array_id,bs='re')",
            "s(array_id,sin1,bs='re')", "s(array_id,cos1,bs='re')",
            "s(array_id,sin2,bs='re')", "s(array_id,cos2,bs='re')", sep=" + ")
marg_s <- if (BASIS=="ds") "bs='ds',k=60,m=c(1,0.5)" else sprintf("bs='gp',k=60,m=c(3,%g)", GPR)
marg_ti <- if (BASIS=="ds") "bs=c('cc','ds'),d=c(1,2),k=c(10,60),m=list(NA,c(1,0.5))" else
                            sprintf("bs=c('cc','gp'),d=c(1,2),k=c(10,60),m=list(NA,c(3,%g))", GPR)
ftxt <- sprintf("count ~ s(bin_center_h,bs='cc',k=10) + s(X,Y,%s) + ti(bin_center_h,X,Y,%s) + %s%s + offset(log_effort_aso)",
                marg_s, marg_ti, RE, if (is.null(EXTRA)) "" else paste0(" + ", EXTRA))
form <- as.formula(ftxt)

## ---- OFFSET ASSERTIONS (before fitting)
stopifnot(grepl("offset\\(log_effort_aso\\)", paste(deparse(form), collapse=" ")))
stopifnot(all(is.finite(d$log_effort_aso)), all(d$effort_h_aso > 0))

t0 <- Sys.time()
m <- bam(form, family=nb(), data=d, discrete=TRUE, method="fREML", nthreads=10,
         knots=list(bin_center_h=c(0, 24)))
mins <- as.numeric(difftime(Sys.time(), t0, units="mins"))
## ---- OFFSET ASSERTION (offset live in the fitted object)
stopifnot(!is.null(m$offset), max(abs(m$offset - d$log_effort_aso)) < 1e-9)

re_lab <- function(m) { L <- sapply(m$smooth, function(s) s$label); L[grepl("array_id", L)] }
edf <- pen.edf(m); grp <- tapply(edf, sub("\\.[0-9]+$","",names(edf)), sum)
pr <- residuals(m, "pearson")
disp <- sum(pr^2)/df.residual(m)
mu <- fitted(m); th <- m$family$getTheta(TRUE)
zero_ratio <- sum(d$count==0)/sum((th/(th+mu))^th)

## cyclic-join check on the population diel main effect
gj <- data.table(bin_center_h=c(1e-4, 24-1e-4), X=median(d$X), Y=median(d$Y),
                 sin1=0, cos1=1, sin2=0, cos2=1, log_effort_aso=0,
                 array_id=d$array_id[1], dep_key=d$dep_key[1])
if (!is.null(EXTRA)) gj[, ecoregion_l1 := d$ecoregion_l1[1]]
tj <- predict(m, gj, type="terms")[, "s(bin_center_h)"]
join_gap <- abs(diff(tj))

## population diel curve vs raw, on the real data rows
lp_pop <- predict(m, d, type="link", exclude=re_lab(m))
agg <- data.table(h=d$bin_center_h, r=exp(lp_pop - d$log_effort_aso))[, .(rate=mean(r)), by=h][order(h)]
raw <- d[, .(rate=sum(count)/sum(effort_h_aso)), by=bin_center_h][order(bin_center_h)]
pn <- function(rate,h){n<-h<SUNRISE|h>SUNSET; 100*sum(rate[n])/sum(rate)}
## mean of array harmonic REs (identifiability check)
cf <- coef(m); remeans <- sapply(c("s(sin1,array_id)","s(cos1,array_id)"), function(tg) {
  i <- which(sapply(m$smooth, function(s) s$label)==tg)
  mean(cf[m$smooth[[i]]$first.para:m$smooth[[i]]$last.para]) })

## ---- grid predictions (population level) + array-level residuals for Moran's I
GRID <- fread("pred_grid_25km.csv")
g <- CJ(i=seq_len(nrow(GRID)), bin_center_h=H)
g[, `:=`(X=GRID$X[i], Y=GRID$Y[i])]
g[, `:=`(sin1=sin(2*pi*bin_center_h/24), cos1=cos(2*pi*bin_center_h/24),
         sin2=sin(4*pi*bin_center_h/24), cos2=cos(4*pi*bin_center_h/24),
         log_effort_aso=0, array_id=d$array_id[1], dep_key=d$dep_key[1])]
if (!is.null(EXTRA)) g[, ecoregion_l1 := d$ecoregion_l1[1]]
LP <- matrix(as.numeric(predict(m, g, type="link", exclude=re_lab(m))), nrow=length(H))
## SE of log rate at each grid cell, averaged over bins (uncertainty surface input)
SE <- matrix(as.numeric(predict(m, g, type="link", exclude=re_lab(m), se.fit=TRUE)$se.fit), nrow=length(H))

arr_res <- data.table(array_id=as.character(d$array_id), X=d$X, Y=d$Y,
                      dres=residuals(m, "deviance"), pres=pr)[
  , .(mean_dres=mean(dres), mean_pres=mean(pres), X=mean(X), Y=mean(Y), n=.N), by=array_id]

out <- list(sp=SP, basis=BASIS, gp_range=GPR, extra=EXTRA, form=gsub("\\s+"," ",ftxt),
            mins=mins, n=nrow(d), ndep=nlevels(d$dep_key), narray=nlevels(d$array_id),
            dev_expl=summary(m)$dev.expl, theta=th, dispersion=disp,
            deviance_disp=sum(residuals(m,"deviance")^2)/df.residual(m),
            converged=isTRUE(m$converged), reml=m$gcv.ubre, aic=AIC(m),
            edf_groups=round(grp,2), edf_total=sum(m$edf), n_coef=length(coef(m)),
            zero_ratio=zero_ratio, join_gap=join_gap, re_means=remeans,
            pop_pct=pn(agg$rate,agg$h), raw_pct=pn(raw$rate,raw$bin_center_h),
            pop_amp=diff(range(log(agg$rate))), raw_amp=diff(range(log(raw$rate))),
            pop_curve=agg, raw_curve=raw, LP=LP, SE=SE, arr_res=arr_res,
            k_check=k.check(m))
saveRDS(out, sprintf("out_%s.rds", TAG))

## slim model retaining only what predict() needs (coefficients, Vp, smooth setup)
slim <- function(mm) {
  for (nm in c("model","residuals","fitted.values","linear.predictors","weights",
               "prior.weights","y","offset","working.weights","z","Xcentre",
               "qrx","R","F","hat","std.rsd")) mm[[nm]] <- NULL
  for (a in c("formula","pterms","terms")) attr(mm[[a]], ".Environment") <- globalenv()
  mm
}
ms <- slim(m)
lp_chk <- as.numeric(predict(ms, g[1:96], type="link", exclude=re_lab(m)))
stopifnot(max(abs(lp_chk - as.numeric(LP)[1:96])) < 1e-9)
saveRDS(ms, sprintf("slim2_%s.rds", TAG))
cat(sprintf("%-22s %-3s%s | %5.2f min dev=%.4f theta=%.3f disp=%.3f join=%.2e | POP pct=%5.1f (raw %5.1f) amp=%.2f (raw %.2f) | REmean cos1=%+.4f\n",
    SP, BASIS, if (!is.na(GPR)) paste0(GPR) else "", mins, out$dev_expl, th, disp,
    join_gap, out$pop_pct, out$raw_pct, out$pop_amp, out$raw_amp, remeans[2]))
