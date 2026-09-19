#!/usr/bin/env python3
"""
submit_chains.py -- the only sanctioned way to submit single-shot iSDM chains.

WHY
Chains were requested at 110G against a measured all-time peak of 32.4G. Hazel
bills memory at 1 unit per 4G, so a 110G chain counts as 28 units of usage, 27
of them memory, and our fair-share factor (0.025 on 9/19) -- not the QOS CPU cap
and not cluster load -- is what made chains wait 30-40 h. This script derives
the memory and wall request from measured runs so the request can't drift back.

WHAT IT DOES (per run dir)
  1. Reads cell50/nsite from the bundle (cached in <run_dir>/bundle_meta.json,
     keyed on the bundle's size+mtime; recomputed with Rscript if stale).
  2. Size class by cell50 -> memory request (small 16G, large 48G).
     REFUSES if the request is > MAX_RATIO x the largest measured peak in that
     class, or < MIN_RATIO x it (OOM protection), or if the class has no
     measured peak, or if cell50 extrapolates > EXTRAP beyond the class's
     largest measured cell50.
  3. Wall = ceil(WALL_MULT x predicted), predicted = cell50 x (hours/cell50)
     from measured completed chains of the same model kind -- same species if
     one exists, else the slowest species of that kind. QOS normal if the wall
     fits its 96 h cap, else long (240 h cap here; refuses beyond it).
  4. Harvests any finished chains into the calibration table first, so every
     submission uses the latest measurements.
  5. Appends one row per chain to chain_submissions.csv.

USAGE
  submit_chains.py plan   <run_dir> [--chains 1,2,3]      # dry run, prints the request
  submit_chains.py submit <run_dir> [--chains 1,2,3]
  submit_chains.py register <run_dir> <chain> <jobid> <req_mem_gb> <req_wall_h> <qos>
                                                           # log a job not submitted here
  submit_chains.py harvest                                 # pull finished chains into calibration
"""
import csv, json, math, os, re, subprocess, sys, datetime

PROJ = "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
HPC = os.path.join(PROJ, "HPC")
HERE = os.path.dirname(os.path.abspath(__file__))
CALIB = os.path.join(HERE, "chain_calibration.csv")
SUBS = os.path.join(HERE, "chain_submissions.csv")
RSCRIPT = os.path.join(HPC, "conda_envs/nimble_env/bin/Rscript")
META_R = os.path.join(HERE, "bundle_meta.R")

# Size classes: (name, request_gb). A run dir is "small" if its cell50 is within
# EXTRAP of the largest measured small-class cell50, else "large".
REQ_GB = {"small": 16, "large": 48}
SMALL_MAX_MEASURED_CELL50_FALLBACK = 526   # moose union; used only if calib has no small rows
MAX_RATIO = 2.0     # request must be <= 2x class peak
MIN_RATIO = 1.3     # request must be >= 1.3x class peak
EXTRAP = 1.15       # refuse if cell50 > 1.15x the largest measured cell50 in class
WALL_MULT = 2.0
QOS_CAP_H = {"normal": 96, "long": 240}
MIN_WALL_H = 4

CALIB_COLS = ["run_dir", "species", "kind", "chain", "jobid", "cell50", "nsite", "req_mem_gb",
              "peak_gb", "wall_h", "state", "source"]
SUB_COLS = ["submitted_at", "run_dir", "species", "kind", "chain", "jobid", "cell50", "nsite", "size_class",
            "req_mem_gb", "class_peak_gb", "req_wall_h", "pred_wall_h", "rate_h_per_cell50", "rate_source",
            "qos", "harvested"]


def die(msg):
    sys.stderr.write("REFUSED: " + msg + "\n")
    sys.exit(3)


def read_csv(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f))


def write_csv(path, cols, rows):
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    os.replace(tmp, path)


def kind_of(run_dir):
    if "camcount" in run_dir:
        return "camcount"
    if run_dir.endswith("_ecoregion") or "_ecoregion_" in run_dir:
        return "ecoregion"
    if "national_scalar" in run_dir:
        return "national_scalar"
    die("cannot tell model kind from run dir name: " + run_dir)


def species_of(run_dir):
    m = re.match(r"^(.*?)_v\d", run_dir)
    if not m:
        die("cannot parse species from run dir name: " + run_dir)
    return m.group(1)


def runner_of(run_dir):
    d = os.path.join(HPC, run_dir)
    for name in ("run_single_camcount_sbatch.sh", "run_single_sbatch.sh"):
        if os.path.exists(os.path.join(d, name)):
            return os.path.join(d, name)
    die("no run_single*_sbatch.sh in " + d)


def bundle_meta(run_dir):
    d = os.path.join(HPC, run_dir)
    rds = os.path.join(d, "input_data_%s.RDS" % run_dir)
    if not os.path.exists(rds):
        die("bundle not found: " + rds)
    st = os.stat(rds)
    key = {"size": st.st_size, "mtime": int(st.st_mtime)}
    cache = os.path.join(d, "bundle_meta.json")
    if os.path.exists(cache):
        with open(cache) as f:
            m = json.load(f)
        if m.get("size") == key["size"] and m.get("mtime") == key["mtime"]:
            return m
    out = subprocess.run([RSCRIPT, META_R, rds], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         universal_newlines=True)
    if out.returncode != 0:
        die("bundle_meta.R failed:\n" + out.stderr)
    m = json.loads(out.stdout.strip().splitlines()[-1])
    m.update(key)
    with open(cache, "w") as f:
        json.dump(m, f)
    return m


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def size_class(cell50, calib):
    small = [fnum(r["cell50"]) for r in calib if fnum(r["cell50"]) and fnum(r["cell50"]) < 1000]
    small_max = max(small) if small else SMALL_MAX_MEASURED_CELL50_FALLBACK
    return ("small" if cell50 <= EXTRAP * small_max else "large"), small_max


def class_rows(cls, small_max, calib):
    out = []
    for r in calib:
        c = fnum(r["cell50"])
        if c is None or fnum(r["peak_gb"]) is None:
            continue
        if (cls == "small") == (c <= EXTRAP * small_max):
            out.append(r)
    return out


def plan(run_dir, chains):
    calib = read_csv(CALIB)
    meta = bundle_meta(run_dir)
    cell50, nsite = int(meta["ncell50"]), int(meta["nsite"])
    kind, sp = kind_of(run_dir), species_of(run_dir)

    cls, small_max = size_class(cell50, calib)
    rows = class_rows(cls, small_max, calib)
    if not rows:
        die("no measured peak in size class %s" % cls)
    peak = max(fnum(r["peak_gb"]) for r in rows)
    max_c50 = max(fnum(r["cell50"]) for r in rows)
    if cell50 > EXTRAP * max_c50:
        die("cell50 %d is > %.2fx the largest measured cell50 in class %s (%d)" % (cell50, EXTRAP, cls, max_c50))
    req = REQ_GB[cls]
    if req > MAX_RATIO * peak:
        die("%dG is > %.1fx the class %s peak of %.1fG" % (req, MAX_RATIO, cls, peak))
    if req < MIN_RATIO * peak:
        die("%dG is < %.1fx the class %s peak of %.1fG (OOM risk)" % (req, MIN_RATIO, cls, peak))

    timed = [r for r in calib if r["kind"] == kind and r["state"] == "COMPLETED"
             and fnum(r["wall_h"]) and fnum(r["cell50"])]
    same = [r for r in timed if r["species"] == sp]
    pool, src = (same, "same species") if same else (timed, "slowest species of kind")
    if not pool:
        die("no completed %s chain with a measured wall" % kind)
    rate = max(fnum(r["wall_h"]) / fnum(r["cell50"]) for r in pool)
    pred = cell50 * rate
    wall = max(MIN_WALL_H, int(math.ceil(WALL_MULT * pred)))
    qos = "normal" if wall <= QOS_CAP_H["normal"] else "long"
    if wall > QOS_CAP_H["long"]:
        die("wall %dh exceeds the long QOS cap" % wall)
    return dict(run_dir=run_dir, species=sp, kind=kind, cell50=cell50, nsite=nsite, size_class=cls,
                req_mem_gb=req, class_peak_gb=round(peak, 1), req_wall_h=wall, pred_wall_h=round(pred, 1),
                rate_h_per_cell50=round(rate, 5), rate_source=src, qos=qos, chains=chains)


def show(p):
    print("run dir     %s  (%s, %s)" % (p["run_dir"], p["species"], p["kind"]))
    print("bundle      cell50 %d  nsite %d" % (p["cell50"], p["nsite"]))
    print("memory      %dG  (class %s, measured peak %.1fG -> %.2fx)" % (
        p["req_mem_gb"], p["size_class"], p["class_peak_gb"], p["req_mem_gb"] / p["class_peak_gb"]))
    print("wall        %dh = %.1fx predicted %.1fh (%.5f h/cell50, %s)" % (
        p["req_wall_h"], WALL_MULT, p["pred_wall_h"], p["rate_h_per_cell50"], p["rate_source"]))
    print("qos         %s   chains %s" % (p["qos"], ",".join(map(str, p["chains"]))))


def submit(p):
    runner = runner_of(p["run_dir"])
    d = os.path.join(HPC, p["run_dir"])
    if os.path.exists(os.path.join(d, "BUNDLE_REVIEW_REQUIRED")):
        die("BUNDLE_REVIEW_REQUIRED present in " + d)
    subs = read_csv(SUBS)
    for c in p["chains"]:
        if os.path.exists(os.path.join(d, "chain_%s_%d_single.RDS" % (p["run_dir"], c))):
            die("chain %d already has a checkpoint in %s" % (c, d))
        cmd = ["sbatch", "--parsable",
               "--job-name=%s_single_c%d" % (p["run_dir"], c),
               "--mem=%dG" % p["req_mem_gb"], "--time=%d:00:00" % p["req_wall_h"], "--qos=" + p["qos"],
               "--output=%s/slurm_single_c%d_%%j.log" % (d, c),
               "--export=ALL,SPECIES=%s,CHAIN_ID=%d,SUBMIT_GUARD=1" % (p["run_dir"], c), runner]
        out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if out.returncode != 0:
            sys.stderr.write("sbatch failed for chain %d: %s\n" % (c, out.stderr.strip()))
            write_csv(SUBS, SUB_COLS, subs)
            sys.exit(4)
        jid = out.stdout.strip().split(";")[0]
        print("chain %d -> job %s" % (c, jid))
        row = {k: p.get(k, "") for k in SUB_COLS}
        row.update(submitted_at=datetime.datetime.now().isoformat(timespec="seconds"), chain=c, jobid=jid,
                   harvested="no")
        subs.append(row)
        write_csv(SUBS, SUB_COLS, subs)


def register(run_dir, chain, jobid, mem, wall, qos):
    meta = bundle_meta(run_dir)
    subs = read_csv(SUBS)
    if any(r["jobid"] == jobid for r in subs):
        print("job %s already registered" % jobid)
        return
    subs.append(dict(submitted_at="registered " + datetime.datetime.now().isoformat(timespec="seconds"),
                     run_dir=run_dir, species=species_of(run_dir), kind=kind_of(run_dir), chain=chain,
                     jobid=jobid, cell50=meta["ncell50"], nsite=meta["nsite"], req_mem_gb=mem,
                     req_wall_h=wall, qos=qos, harvested="no"))
    write_csv(SUBS, SUB_COLS, subs)
    print("registered job %s" % jobid)


def to_gb(s):
    s = s.strip()
    if not s:
        return None
    mult = {"K": 1 / 1024 ** 2, "M": 1 / 1024, "G": 1, "T": 1024}
    return float(s[:-1]) * mult[s[-1]] if s[-1] in mult else float(s) / 1024 ** 3


def to_h(s):
    d = 0
    if "-" in s:
        d, s = s.split("-")
        d = int(d)
    p = [int(x) for x in s.split(":")]
    while len(p) < 3:
        p.insert(0, 0)
    return d * 24 + p[0] + p[1] / 60 + p[2] / 3600


def harvest(quiet=False):
    subs, calib = read_csv(SUBS), read_csv(CALIB)
    pending = [r for r in subs if r["harvested"] == "no"]
    if not pending:
        return
    ids = ",".join(r["jobid"] for r in pending)
    out = subprocess.run(["sacct", "-n", "-P", "-j", ids, "-o", "JobID,State,Elapsed,MaxRSS"],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if out.returncode != 0:
        sys.stderr.write("harvest: sacct failed, skipping: %s\n" % out.stderr.strip())
        return
    state, elapsed, rss = {}, {}, {}
    for line in out.stdout.splitlines():
        jid, st, el, mr = (line.split("|") + ["", "", "", ""])[:4]
        base = jid.split(".")[0]
        if "." not in jid:
            state[base], elapsed[base] = st.split()[0], el
        elif jid.endswith(".batch"):
            rss[base] = to_gb(mr)
    n = 0
    known = {r["jobid"] for r in calib}
    for r in pending:
        j = r["jobid"]
        st = state.get(j, "")
        if st in ("", "PENDING", "RUNNING", "REQUEUED", "SUSPENDED"):
            continue
        if j in known:          # already in calibration (e.g. seeded); just mark it
            r["harvested"] = "yes"
            n += 1
            continue
        calib.append(dict(run_dir=r["run_dir"], species=r["species"], kind=r["kind"], chain=r["chain"], jobid=j,
                          cell50=r["cell50"], nsite=r["nsite"], req_mem_gb=r["req_mem_gb"],
                          peak_gb=round(rss[j], 2) if rss.get(j) else "",
                          wall_h=round(to_h(elapsed[j]), 2) if st == "COMPLETED" else "",
                          state=st, source="harvest"))
        r["harvested"] = "yes"
        n += 1
        if not quiet:
            print("harvested %s %s c%s: %s peak %s G wall %s h" % (j, r["run_dir"], r["chain"], st,
                  calib[-1]["peak_gb"], calib[-1]["wall_h"]))
    if n:
        write_csv(CALIB, CALIB_COLS, calib)
        write_csv(SUBS, SUB_COLS, subs)


def parse_chains(args):
    if "--chains" in args:
        return [int(x) for x in args[args.index("--chains") + 1].split(",")]
    return [1, 2, 3]


if __name__ == "__main__":
    a = sys.argv[1:]
    if not a:
        print(__doc__)
        sys.exit(1)
    if a[0] == "harvest":
        harvest()
    elif a[0] in ("plan", "submit"):
        harvest(quiet=True)
        p = plan(a[1], parse_chains(a))
        show(p)
        if a[0] == "submit":
            submit(p)
    elif a[0] == "register":
        register(*a[1:7])
    else:
        print(__doc__)
        sys.exit(1)
