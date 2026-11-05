#!/usr/bin/env bash
#SBATCH --job-name=win_smoke
#SBATCH --partition=compute
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --output=/home/rwkays/isdm/fleet_v2u/windowed/logs/win_smoke_%x_%j.log
# Smoke test for a windowed run dir (2026-10-17). NOT a reported chain: runs the run dir's own
# HPC_run_model_single.R for NBURNIN=100 + NITER=300 in a throwaway copy HPC/<RD>_smoke/ (runner files
# copied, bundle symlinked), so the real chain outputs in HPC/<RD>/ stay untouched and the guard still
# owns every reported chain. Checks: model builds and compiles on the windowed constants, the MCMC runs,
# year_beta / total_var_beta are sampled and finite. Chain id 12 (exercises the 1:12 chain cap).
# Usage: sbatch --job-name=sm_<RD> --mem=<M> --time=<T> --export=ALL,RD=<run_dir> run_window_smoke.sh
set -uo pipefail
: "${RD:?set RD}"
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
NE="$PROJ/HPC/conda_envs/nimble_env"; export PATH="$NE/bin:$PATH"
S=$PROJ/HPC/$RD; T=$PROJ/HPC/${RD}_smoke
[ -f "$S/HPC_run_model_single.R" ] || { echo "FAIL: $S not populated"; exit 1; }
[ -e "$T" ] && { echo "FAIL: $T exists"; exit 1; }
[ -e "$S/BUNDLE_REVIEW_REQUIRED" ] && { echo "FAIL: $S carries BUNDLE_REVIEW_REQUIRED"; exit 1; }
mkdir -p $T
cp $S/HPC_run_model_chunks_chain1.R $S/HPC_run_model_single.R $S/pipeline_helper.R $S/integration_helper_v2.R $T/
ln -s $S/input_data_${RD}.RDS $T/input_data_${RD}_smoke.RDS
echo "host $(hostname) start $(date) RD $RD"
cd $T && SPECIES=${RD}_smoke CHAIN_ID=12 NITER=300 NBURNIN=100 "$NE/bin/Rscript" HPC_run_model_single.R; RC=$?
[ $RC = 0 ] && "$NE/bin/Rscript" -e "x <- readRDS('chain_${RD}_smoke_12_single.RDS'); s <- x\$samples
  p <- intersect(c('year_beta','total_var_beta','year_var','theta1'), colnames(s))
  print(t(apply(s[, p, drop = FALSE], 2, function(v) c(mean = mean(v), sd = sd(v), n_unique = length(unique(v))))))
  stopifnot('year_beta' %in% p, all(is.finite(s[, p])), length(unique(s[, 'year_beta'])) > 1)
  cat('build_sec', round(x\$build_sec), ' sec_per_iter', signif(x\$mcmc_sec / 400, 3), '\nSMOKE OK $RD\n')"
RC2=$?; echo "exit $RC/$RC2 end $(date)"; [ $RC = 0 ] && [ $RC2 = 0 ]
