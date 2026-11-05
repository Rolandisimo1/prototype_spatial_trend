#!/usr/bin/env bash
#SBATCH --job-name=win_feeder
#SBATCH --partition=compute
#SBATCH --qos=long
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/home/rwkays/isdm/fleet_v2u/windowed/logs/feeder_%j.log
# QOS long (2026-10-19): a normal job is refused once ALL compute jobs reach 1000, which would kill the
# re-queue; long has its own count and the feeder keeps it under CAP 980.
# One feeder pass (window_feeder.py), then re-queue itself to start >= 3 h later while run dirs remain.
# Stop it with: touch ~/isdm/fleet_v2u/windowed/FEEDER_STOP  (checked before each pass and before re-queueing)
set -uo pipefail
WD=/home/rwkays/isdm/fleet_v2u/windowed; Q=$WD/window_queue.txt
[ -e $WD/FEEDER_STOP ] && { echo "FEEDER_STOP present -- exiting"; exit 0; }
python3 $WD/window_feeder.py $Q
LEFT=$(cat $Q.left 2>/dev/null || echo 1)
if [ "$LEFT" -gt 0 ] && [ ! -e $WD/FEEDER_STOP ]; then
  sbatch --begin=now+3hours $WD/window_feeder_sbatch.sh
else echo "feeder finished or stopped: $LEFT run dirs still waiting"; fi
