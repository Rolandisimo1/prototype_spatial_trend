"""Audit species names against the Mammal Diversity Database (MDD).

Why this exists, and what the existing crosswalk could not do.

`fleet_union_prep/build_name_crosswalk.R` reconciles names ACROSS SOURCES: one
row per fleet species, one name per source (inat_name, iucn_name, wi_name,
mdd_binomial). That is the right object for "what does each database call this
species" and it caught the armadillo and fisher cases.

It cannot catch a single source using TWO names over time, because there is only
one slot per source. Wildlife Insights filed the peccary as `Pecari tajacu`
through 2023 and `Dicotyles tajacu` from 2024, and a one-slot schema cannot
represent that. A pipeline keyed on the recorded name silently lost 47.9% of the
record, all of it at the recent end, which read as a -54 %/yr decline.

So this module audits WITHIN a source: resolve every distinct name a source uses
to its MDD species id, and group by id rather than by string. Two names with one
id are a rename to be unioned; two names with two ids are two species and must
NOT be unioned.

Pins to an MDD version and records it, per the mdd-taxonomy skill. Comparison is
on id, never on name strings.
"""
import os

import numpy as np
import pandas as pd

MDD_DIR = "/Users/rwkays/Documents/iNat/Global_inat/Datasets/otherdata/MDD"
MDD_VERSION = "v2.5"

# Human decisions that name resolution CANNOT make. Two MDD-accepted names are
# two species by definition, so merging them is a source-data correction, not a
# lookup -- it needs an authority, a basis and a date, recorded here rather than
# hardcoded. Same pattern as manual_range_overrides.csv.
DECISIONS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                         "fleet_union_prep", "manual_taxonomy_decisions.csv")


def load_mdd(mdd_dir=MDD_DIR, version=MDD_VERSION):
    """Accepted-name -> id map and synonym -> accepted-name map, both from MDD tables."""
    spf = [f for f in os.listdir(mdd_dir)
           if f.startswith(f"MDD_{version}") and f.endswith(".csv")]
    synf = [f for f in os.listdir(mdd_dir)
            if f.startswith("Species_Syn_Current") and version in f]
    if not spf or not synf:
        raise FileNotFoundError(f"MDD {version} species/synonym CSVs not found in {mdd_dir}")
    sp = pd.read_csv(os.path.join(mdd_dir, spf[0]), low_memory=False)
    syn = pd.read_csv(os.path.join(mdd_dir, synf[0]), low_memory=False)
    for col, src in [("sciName", sp), ("MDD_original_combination", syn)]:
        if col not in src.columns:
            raise KeyError(f"expected column {col!r} absent; MDD schema may have changed")
    sp["sciName_n"] = sp.sciName.str.replace("_", " ").str.strip()
    # Genus-specific match key: the full original combination, never a bare
    # epithet -- epithets are reused across unrelated genera (mdd-taxonomy trap).
    syn["orig_n"] = syn.MDD_original_combination.astype(str).str.replace("_", " ").str.strip()
    return (dict(zip(sp.sciName_n, sp.id)),
            syn.dropna(subset=["MDD_species"]).set_index("orig_n").MDD_species.to_dict(),
            version)


def resolve(names, mdd_dir=MDD_DIR, version=MDD_VERSION):
    """Resolve names to MDD ids. Unresolved stays unresolved -- never invented."""
    acc_id, syn_map, ver = load_mdd(mdd_dir, version)
    out = []
    for n in pd.unique(pd.Series(list(names)).dropna().astype(str).str.strip()):
        if n in acc_id:
            out.append({"verbatim": n, "resolved": n, "mdd_id": acc_id[n],
                        "match_type": "accepted_exact"})
            continue
        s = str(syn_map.get(n, "")).replace("_", " ").strip()
        if s and s in acc_id:
            out.append({"verbatim": n, "resolved": s, "mdd_id": acc_id[s],
                        "match_type": "synonym"})
        else:
            out.append({"verbatim": n, "resolved": np.nan, "mdd_id": np.nan,
                        "match_type": "unresolved"})
    d = pd.DataFrame(out)
    d["mdd_version"] = ver
    return d


def audit_within_source(df, name_col, group_col, mdd_dir=MDD_DIR,
                        version=MDD_VERSION, decisions=None):
    """Find group labels (e.g. common names) carrying more than one scientific name.

    Returns one row per such group, with `n_mdd_ids` deciding the treatment:
      1  -> a rename; union the names, they are one species.
      >1 -> distinct species sharing a label; unioning them would merge taxa.
            Needs a documented decision, not a lookup. If one is recorded in
            manual_taxonomy_decisions.csv it is reported in the `decision`
            column -- the flag is never silently cleared, so a reader always
            sees that a human made the call and who.
    """
    res = resolve(df[name_col].dropna().unique(), mdd_dir, version).set_index("verbatim")
    d = df[[group_col, name_col]].dropna().drop_duplicates().copy()
    d["mdd_id"] = d[name_col].astype(str).str.strip().map(res.mdd_id)
    g = d.groupby(group_col).agg(
        names=(name_col, lambda s: " | ".join(sorted(set(map(str, s))))),
        n_names=(name_col, "nunique"), n_mdd_ids=("mdd_id", "nunique"),
        n_unresolved=("mdd_id", lambda s: int(s.isna().sum())))
    g = g[(g.n_names > 1) | (g.n_unresolved > 0)].reset_index()
    g["treatment"] = np.where(g.n_unresolved > 0, "UNRESOLVED - resolve before use",
                       np.where(g.n_mdd_ids == 1, "rename - union the names",
                                "DISTINCT SPECIES - needs a recorded decision"))
    g["decision"] = ""
    if os.path.exists(decisions or DECISIONS):
        dec = pd.read_csv(decisions or DECISIONS)
        for _, r in dec.iterrows():
            m = g[group_col].astype(str).str.lower() == str(r.common_name).lower()
            if not m.any():
                continue
            g.loc[m, "decision"] = (f"{r.from_name} -> {r.to_name} "
                                    f"(id {r.to_mdd_id}); {r.authority}, {r.date}")
            g.loc[m & g.treatment.str.startswith("DISTINCT"), "treatment"] = (
                "DISTINCT SPECIES - resolved by recorded decision")
    g["mdd_version"] = version
    return g


def audit_crosswalk(crosswalk, mdd_col="mdd_binomial", mdd_dir=MDD_DIR, version=MDD_VERSION):
    """Check a crosswalk's own MDD column actually holds MDD-accepted names.

    The mdd-taxonomy checklist requires every match to come from a table lookup.
    A row whose mdd_binomial is not in the accepted list means that lookup did
    not happen, or was not verified, for that row. One line, catches it.
    """
    acc_id, syn_map, ver = load_mdd(mdd_dir, version)
    d = crosswalk.copy()
    d["_n"] = d[mdd_col].astype(str).str.strip()
    d["is_mdd_accepted"] = d._n.isin(acc_id)
    d["should_be"] = d._n.map(lambda x: str(syn_map.get(x, "")).replace("_", " ").strip())
    d["mdd_version"] = ver
    return d.loc[~d.is_mdd_accepted,
                 [c for c in d.columns if c != "_n"]].reset_index(drop=True)
