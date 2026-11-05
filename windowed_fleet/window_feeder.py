#!/usr/bin/env python3
"""window_feeder.py (2026-10-17) -- staged launch of the windowed refits under the per-user QOS submit caps.

Hazel caps a user at 1000 SUBMITTED jobs per QOS (normal and long each; every array task counts). The windowed
refits are 136 fits x 12 chains = 1632 jobs, alongside Tier 1 (540 array tasks) and the ybfix sweep, so they
cannot all be queued at once. Each pass: read the queue file (one run dir per line, in launch order), and for
each run dir not yet submitted ask the guard for its plan; if that plan's QOS has room for all 12 chains under
CAP (submitted jobs counted with squeue -r), submit it THROUGH THE GUARD (submit_chains.py submit ... --chains
1..12). Every reported chain still goes through submit_chains.py; this script only decides WHEN.
State: <queue>.done (run dir, qos|PARTIAL, job ids; a partial run dir is topped up with its missing chains on a\nlater pass), <queue>.left, feeder.log. Never resubmits a chain recorded in .done. 2026-10-19: a MaxSubmit refusal\nends the pass (the real limit is a per-user total across QOS).
usage: window_feeder.py <queue_file>   (one pass; window_feeder_sbatch.sh re-queues itself every 3 h)
"""
import os, re, subprocess, sys, datetime

GUARD = "/home/rwkays/isdm/fleet_v2u/submit_guard/submit_chains.py"
CAP = 980                        # leave headroom under the 1000 MaxSubmitPU per QOS
CHAINS = ",".join(str(c) for c in range(1, 13))


def log(msg):
    line = "%s %s" % (datetime.datetime.now().isoformat(timespec="seconds"), msg)
    print(line)
    with open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[1])), "feeder.log"), "a") as fh:
        fh.write(line + "\n")


def submitted_by_qos():
    out = subprocess.run(["squeue", "-u", os.environ.get("USER", "rwkays"), "-h", "-r", "-o", "%q"],
                         stdout=subprocess.PIPE, universal_newlines=True, check=True).stdout.split()
    # Partition `compute` carries QOS normal as its PARTITION QOS (OverPartQOS on long), so a normal submission
    # is checked against ALL of the user's compute jobs (normal + long), while long is checked against long only.
    return {"long": out.count("long"), "normal": len(out)}


def read_done(donef):
    """rd -> job ids submitted so far (summed over its .done lines; a PARTIAL line plus a later TOPUP line)."""
    got = {}
    if os.path.exists(donef):
        for l in open(donef):
            f = l.split()
            if len(f) >= 3:
                got.setdefault(f[0], []).extend(f[2].split(","))
    return got


def main(qfile):
    todo = [l.strip() for l in open(qfile) if l.strip() and not l.startswith("#")]
    donef = qfile + ".done"
    got = read_done(donef)
    n = submitted_by_qos()
    log("pass: %d queued, %d complete; submitted jobs by QOS %s" % (
        len(todo), sum(len(got.get(rd, [])) >= 12 for rd in todo), n))
    for rd in todo:
        have = len(got.get(rd, []))
        if have >= 12:
            continue
        # the guard submits chains in order, so a partial run dir holds chains 1..have
        chains = ",".join(str(c) for c in range(have + 1, 13))
        need = 12 - have
        p = subprocess.run(["python3", GUARD, "plan", rd, "--chains", chains], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, universal_newlines=True)
        if p.returncode != 0:
            log("PLAN REFUSED %s: %s" % (rd, p.stderr.strip().splitlines()[-1] if p.stderr.strip() else p.returncode))
            continue
        m = re.search(r"^qos\s+(\S+)", p.stdout, re.M)
        qos = m.group(1)
        if n.get(qos, 0) + need > CAP:
            log("cap: %s holds %d submitted; %s waits" % (qos, n.get(qos, 0), rd))
            continue                                # a later run dir may fit the other QOS
        s = subprocess.run(["python3", GUARD, "submit", rd, "--chains", chains], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, universal_newlines=True)
        jobs = re.findall(r"job (\d+)", s.stdout)
        if jobs:
            with open(donef, "a") as fh:
                fh.write("%s %s %s\n" % (rd, qos if len(jobs) == need else "PARTIAL", ",".join(jobs)))
            n[qos] = n.get(qos, 0) + len(jobs)
            if qos == "long":
                n["normal"] += len(jobs)            # long jobs also fill the partition QOS
        if s.returncode != 0 or len(jobs) != need:
            log("SUBMIT FAILED %s chains %s rc %d jobs %s: %s" % (rd, chains, s.returncode, jobs, s.stderr.strip()[-300:]))
            if "MaxSubmit" in s.stderr:
                # Hazel's binding limit is a per-USER total (~1000 across QOS), not per QOS:
                # nothing else will go in this pass.
                log("submit limit reached; ending pass")
                break
            continue
        log("SUBMITTED %s qos %s chains %s jobs %s-%s" % (rd, qos, chains, jobs[0], jobs[-1]))
    got = read_done(donef)
    left = [rd for rd in todo if len(got.get(rd, [])) < 12]
    log("pass end: %d run dirs still waiting" % len(left))
    with open(qfile + ".left", "w") as fh:
        fh.write("%d\n" % len(left))
    return len(left)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__); sys.exit(1)
    sys.exit(0 if main(sys.argv[1]) >= 0 else 1)
