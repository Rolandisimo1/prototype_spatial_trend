#!/usr/bin/env python3
"""Populate windowed joint run dirs (2026-10-17). usage: build_window_rundirs.py <wtag> <species>:<tag> [...]

Target: HPC/<sp>_<tag>_national_scalar_<wtag>/, which already holds the windowed bundle from
~/isdm/fleet_v2u/windowed/run_window_bundle.sh. Runner files are copied from the FITTED 18-yr run dir
HPC/<sp>_<tag>_national_scalar/ (read-only). One edit, not in the model:
  1. HPC_run_model_chunks_chain1.R: copied verbatim (no token retarget; the single runner takes the
     run dir from SPECIES and reads only the model block + monitors); model block checked identical.
  2. HPC_run_model_single.R: `chain_id %in% 1:3` -> `chain_id %in% 1:12` (rwkays 2026-10-17: 12 chains per fit).
     Seeds stay set.seed(1000 + chain_id); NITER/NBURNIN defaults (50000/5000) unchanged.
The camera-only arm is then built from this dir with
  SRC_TAG=<tag> SRC_SUFFIX=_<wtag> build_variant_rundirs.py camonly <sp>
"""
import hashlib, os, shutil, sys
sys.path.insert(0, "/home/rwkays/isdm/fleet_v2u")
from build_variant_rundirs import nimble_block, FILES, HPC

CAP_OLD = "chain_id %in% 1:3)"
CAP_NEW = "chain_id %in% 1:12)"


def build(wtag, sp, tag):
    src_tok = "%s_%s_national_scalar" % (sp, tag)
    tok = src_tok + "_" + wtag
    src, dst = os.path.join(HPC, src_tok), os.path.join(HPC, tok)
    if not os.path.exists(os.path.join(dst, "input_data_%s.RDS" % tok)):
        return "FAIL %s: windowed bundle missing" % tok
    if any(f.startswith("chain_") for f in os.listdir(dst)):
        return "FAIL %s: holds chain outputs" % tok
    if os.path.exists(os.path.join(dst, "HPC_run_model_single.R")):
        return "SKIP %s: runner already present" % tok
    for f in FILES:
        shutil.copy2(os.path.join(src, f), os.path.join(dst, f))
    # The chunk script is copied VERBATIM (2026-10-18). HPC_run_model_single.R reads only its
    # model_code block and addMonitors list; its `species <-` line is dead text (the fitted v2u
    # dirs still carry the template's "bobcat_v2b_national_scalar"), so it is not retargeted.
    ck = os.path.join(dst, "HPC_run_model_chunks_chain1.R")
    s = open(ck).read()
    so = open(os.path.join(src, "HPC_run_model_chunks_chain1.R")).read()
    a, b = nimble_block(so, "model_code_national_scalar"), nimble_block(s, "model_code_national_scalar")
    if so[a[0]:a[1]] != s[b[0]:b[1]]:
        return "FAIL %s: model block differs from source" % tok
    sg = os.path.join(dst, "HPC_run_model_single.R")
    t = open(sg).read()
    if t.count(CAP_OLD) != 1:
        return "FAIL %s: chain cap anchor not found exactly once" % tok
    open(sg, "w").write(t.replace(CAP_OLD, CAP_NEW))
    with open(os.path.join(dst, "WINDOW_README.txt"), "w") as fh:
        fh.write("windowed run dir %s populated %s from %s (read-only source)\n%s\n" % (
            tok, __import__("datetime").datetime.now().isoformat(timespec="seconds"), src_tok, __doc__))
    md = hashlib.md5(s[b[0]:b[1]].encode()).hexdigest()[:12]
    return "BUILT %s  model block %s (== %s), chunk script verbatim, chain cap 1:12" % (tok, md, src_tok)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    for st in sys.argv[2:]:
        sp, tag = st.split(":")
        print(build(sys.argv[1], sp, tag))
