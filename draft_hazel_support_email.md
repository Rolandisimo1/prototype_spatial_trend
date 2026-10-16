To: help@ncsu.edu
Subject: Hazel Slurm: user rwkays pending jobs get no start estimate despite idle capacity (nodes MIXED+PLANNED)

Hello HPC team,

Since about 06:35 EDT on Sat 10 Oct 2026, none of my jobs on Hazel's compute partition has started, although the
partition has had free capacity throughout. I've ruled out the causes I can see from the user side and would
appreciate a look at the scheduler side.

What I see (checked 18:20-18:45 EDT, 10 Oct):
- 28 jobs pending, account rwkays_cpu, QOS normal, partition compute, reason "Priority". These appear to be the
  only pending jobs on the cluster (`squeue -h -t PD | wc -l` = 28, all mine).
- `squeue -u rwkays --start` reports START_TIME N/A for all 28, including three small jobs
  (1 CPU, 24-48 GB, wall limit lowered to 30-45 min): 1183907, 1183908, 1189006.
- `sinfo -p compute -o "%C"` showed 2440 allocated, 268 idle CPUs. Several nodes have >100 GB free memory, e.g.
  c025n02: CPUAlloc=24/32, FreeMem=168205 MB, State=MIXED+PLANNED.
- `scontrol show reservation` returns "No reservations in the system".
- No GrpTRES/MaxTRES/MaxJobs limits on my association or account. QOS normal: MaxTRESPU cpu=1024, MaxWall 4-00:00:00.
  I have 6 jobs running (6 CPUs).
- Holding my highest-priority pending job (1178027) for five minutes did not let the others start.
- Other users' jobs have kept starting normally during this period.

So the nodes that could take my jobs are marked PLANNED, but no reservation or visible pending job accounts for
that, and the scheduler gives none of my jobs an estimated start. Is something reserving those nodes that isn't
visible to users, such as a maintenance window or a node set held back for the LSF retirement? Or is something
specific to my association keeping my jobs out of scheduling and backfill?

Thank you,
Roland Kays (rwkays)
