## Per-species diagnostic figure: diel curves at contrasting locations, QQ, variogram, ds-vs-gp map
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(viridis)})
SUNRISE <- 6.2202; SUNSET <- 17.9998; H <- seq(0.25,23.75,by=0.5)
tg <- function(s) gsub("[^A-Za-z0-9]+","_",s)
GRID <- fread("pred_grid_25km.csv")
BEST <- readRDS("best_gp.rds")
pctn <- function(LP){r<-exp(LP);p<-sweep(r,2,colSums(r),"/");n<-H<SUNRISE|H>SUNSET;100*colSums(p[n,,drop=FALSE])}

## example locations, chosen as extremes of the fitted surface + a central cell
pick_locs <- function(pn) {
  i <- c(which.min(pn), which.max(pn), which.min(abs(pn - median(pn))))
  data.table(idx=i, label=c("lowest % nocturnal","highest % nocturnal","median"),
             X=GRID$X[i], Y=GRID$Y[i], pct=pn[i])
}

fig_for <- function(sp) {
  o  <- readRDS(sprintf("out_%s_ds.rds", tg(sp)))
  og <- readRDS(sprintf("out_%s_gp_%d.rds", tg(sp), BEST[[sp]]))
  pn <- pctn(o$LP); png <- pctn(og$LP)
  L <- pick_locs(pn)
  cur <- rbindlist(lapply(seq_len(nrow(L)), function(j)
    data.table(h=H, rate=exp(o$LP[, L$idx[j]]), label=sprintf("%s (%.0f%%)", L$label[j], L$pct[j]))))
  p1 <- ggplot(cur, aes(h, rate, colour=label)) +
    annotate("rect", xmin=-Inf, xmax=SUNRISE, ymin=0, ymax=Inf, alpha=.10) +
    annotate("rect", xmin=SUNSET, xmax=Inf, ymin=0, ymax=Inf, alpha=.10) +
    geom_line(linewidth=.8) + scale_x_continuous(breaks=seq(0,24,6), limits=c(0,24)) +
    scale_colour_viridis_d(end=.85) +
    labs(x="sun-anchored hour", y="rate (events / h)", title="a  fitted diel curves", colour=NULL) +
    theme_bw(9) + theme(legend.position="bottom", legend.direction="vertical", legend.key.height=unit(8,"pt"))

  ## QQ of array-level mean deviance residuals
  ar <- o$arr_res[is.finite(mean_dres)]
  qq <- data.table(theo=qnorm(ppoints(nrow(ar))), samp=sort(scale(ar$mean_dres)[,1]))
  p2 <- ggplot(qq, aes(theo, samp)) + geom_abline(slope=1, intercept=0, colour="grey50") +
    geom_point(size=.7, alpha=.7) +
    labs(x="theoretical quantile", y="standardised array residual", title="b  array residual QQ") + theme_bw(9)

  ## variogram of array-level residuals
  v <- readRDS("vario_ds.rds")[[sp]]
  p3 <- ggplot(v, aes(dist, gamma)) + geom_point(aes(size=np), alpha=.7) + geom_line(alpha=.5) +
    scale_size_area(max_size=3, guide="none") + ylim(0, NA) +
    labs(x="separation (km)", y="semivariance", title="c  variogram, array residuals") + theme_bw(9)

  ## ds vs gp scatter of % nocturnal
  dd <- data.table(ds=pn, gp=png, nn=GRID$nn_km)
  p4 <- ggplot(dd, aes(ds, gp, colour=nn)) + geom_abline(slope=1,intercept=0,colour="grey50") +
    geom_point(size=.5, alpha=.5) + scale_colour_viridis_c(name="km to\ncamera", option="magma", end=.9) +
    labs(x="% nocturnal (Duchon)", y=sprintf("%% nocturnal (GP %d km)", BEST[[sp]]),
         title="d  basis comparison") + theme_bw(9) + theme(legend.key.width=unit(6,"pt"))

  ## map of the Duchon surface
  mp <- data.table(X=GRID$X, Y=GRID$Y, pct=pn)
  p5 <- ggplot(mp, aes(X, Y, fill=pct)) + geom_raster() + coord_equal() +
    scale_fill_viridis_c(name="% noct", option="cividis") +
    labs(x=NULL, y=NULL, title="e  fitted % nocturnal (Duchon)") + theme_bw(9) +
    theme(axis.text=element_blank(), axis.ticks=element_blank(), legend.key.width=unit(6,"pt"))

  (p1 | p2 | p3) / (p4 | p5) +
    patchwork::plot_annotation(title=sprintf("%s — space x time-of-day model diagnostics", sp),
      subtitle=sprintf("n=%s rows, %s deployments, %s arrays | dev.expl=%.3f, theta=%.3f, dispersion=%.2f",
        format(o$n,big.mark=","), format(o$ndep,big.mark=","), o$narray, o$dev_expl, o$theta, o$dispersion),
      theme=theme(plot.title=element_text(size=11, face="bold"), plot.subtitle=element_text(size=8)))
}
for (sp in c("White-tailed Deer","Eastern Gray Squirrel","Northern Raccoon","Coyote","American Black Bear")) {
  g <- fig_for(sp)
  ggsave(sprintf("diag_%s.png", tg(sp)), plot=g, width=11, height=7.2, dpi=190)
  cat("wrote diag_", tg(sp), ".png\n", sep="")
}
