## Array-level residual spatial autocorrelation: Moran's I + empirical variogram.
suppressPackageStartupMessages({library(data.table); library(spdep); library(gstat); library(sp)})

moran_of <- function(ar, k=8, nsim=999) {
  ar <- ar[is.finite(mean_dres) & is.finite(X) & is.finite(Y)]
  ## collapse arrays sharing an identical centroid (would give zero-distance neighbours)
  ar <- ar[, .(mean_dres=weighted.mean(mean_dres, n), n=sum(n)), by=.(X, Y)]
  xy <- as.matrix(ar[, .(X, Y)])
  nb <- knn2nb(knearneigh(xy, k=min(k, nrow(ar)-1)), sym=TRUE)
  lw <- nb2listw(nb, style="W")
  mt <- moran.test(ar$mean_dres, lw, randomisation=TRUE)
  mc <- moran.mc(ar$mean_dres, lw, nsim=nsim)
  list(I=as.numeric(mt$estimate["Moran I statistic"]),
       expected=as.numeric(mt$estimate["Expectation"]),
       p_analytic=mt$p.value, p_mc=mc$p.value, n=nrow(ar), k=k)
}

vario_of <- function(ar, cutoff=1500, width=50) {
  ar <- ar[is.finite(mean_dres)][, .(mean_dres=weighted.mean(mean_dres,n)), by=.(X,Y)]
  sp_df <- ar[, .(X, Y, r=mean_dres)]
  coordinates(sp_df) <- ~X+Y
  v <- variogram(r ~ 1, sp_df, cutoff=cutoff, width=width)
  as.data.table(v)[, .(dist, gamma, np)]
}
