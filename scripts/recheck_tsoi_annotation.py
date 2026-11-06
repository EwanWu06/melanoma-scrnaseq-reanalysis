"""Recheck Tsoi-state annotations against the project's own marker table.
See docs/stage4_annotation_recheck.md and decision log entry (j). Read-only on data/processed/."""
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

BASE = Path(__file__).resolve().parents[1] / "data" / "processed"
OUT = Path(__file__).resolve().parents[1] / "results" / "tables" / "recheck"
OUT.mkdir(parents=True, exist_ok=True)
MK = ["MITF", "SOX10", "NGFR", "AXL"]
TABLE = {"Melanocytic": "UUDD", "Transitory": "UUUD", "Neural crest-like": "DUUU", "Undifferentiated": "DDDU"}
lines = []


def log(s=""):
    print(s)
    lines.append(str(s))


def expr_matrix(a):
    idx = [list(a.var_names).index(g) for g in MK]
    X = a.X[:, idx]
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    return pd.DataFrame(X, index=a.obs_names, columns=MK)


def pick_key(a, label_col, keys):
    """The resolution at which the final label is constant within every cluster (fewest clusters wins)."""
    ok = []
    for k in keys:
        if k in a.obs and a.obs.groupby(k, observed=True)[label_col].nunique().max() == 1:
            ok.append((a.obs[k].nunique(), k))
    return sorted(ok)[0][1] if ok else None


def directions(pct, overall, rule):
    d = {}
    for g in MK:
        if rule == 1:
            d[g] = "U" if pct[g] >= 1.25 * overall[g] else ("D" if pct[g] <= 0.75 * overall[g] else "~")
        else:
            d[g] = "U" if pct[g] >= 0.5 else "D"
    return d


def matches(d, genes):
    out = []
    for state, pat in TABLE.items():
        if all(d[g] in ("~", pat[MK.index(g)]) for g in genes):
            out.append(state)
    return out


def analyse(name, fname, label_col, patient_col, keys):
    a = ad.read_h5ad(BASE / fname)
    e = expr_matrix(a)
    log(f"\n==================== {name} ({fname}) n={a.n_obs}")
    log(f"min value of the 4 markers: {e.values.min():.3f} (must be >= 0, i.e. not scaled)")
    key = pick_key(a, label_col, keys)
    log(f"cluster key used for the original labels: {key}")
    pos = e > 0
    overall = pos.mean()
    log("all malignant cells, detection rate: " + ", ".join(f"{g} {overall[g]:.1%}" for g in MK))
    rows = []
    for c, idx in a.obs.groupby(key, observed=True).groups.items():
        pct = pos.loc[idx].mean()
        pats = a.obs.loc[idx, patient_col].astype(str).value_counts()
        d1, d2 = directions(pct, overall, 1), directions(pct, overall, 2)
        rows.append({
            "cluster": c, "n": len(idx), "original_label": a.obs.loc[idx, label_col].iloc[0],
            **{f"{g}%": round(100 * pct[g], 1) for g in MK},
            "rule1": "".join(d1[g] for g in MK), "rule1_match4": "/".join(matches(d1, MK)) or "none",
            "rule1_match3": "/".join(matches(d1, MK[:3])) or "none",
            "rule2": "".join(d2[g] for g in MK), "rule2_match4": "/".join(matches(d2, MK)) or "none",
            "rule2_match3": "/".join(matches(d2, MK[:3])) or "none",
            "top_patient": f"{pats.index[0]} {pats.iloc[0] / len(idx):.0%}",
            "single_patient": bool(pats.iloc[0] / len(idx) >= 0.75),
        })
    t = pd.DataFrame(rows)
    pd.set_option("display.width", 250)
    log(t.to_string(index=False))
    t.to_csv(OUT / f"clusters_{name}.csv", index=False)
    # cell level
    cats = {
        "T-like (NGFR+ SOX10+ MITF+)": pos.NGFR & pos.SOX10 & pos.MITF,
        "NC-like (NGFR+ SOX10+ MITF-)": pos.NGFR & pos.SOX10 & ~pos.MITF,
        "NGFR+ SOX10-": pos.NGFR & ~pos.SOX10,
    }
    cell_rows = []
    for k, m in cats.items():
        pats = a.obs.loc[m.values, patient_col].astype(str).value_counts()
        cell_rows.append({"category": k, "cells": int(m.sum()), "pct_of_malignant": round(100 * m.mean(), 1),
                          "patients_ge3": int((pats >= 3).sum()),
                          "top_patient": f"{pats.index[0]} {pats.iloc[0] / max(m.sum(), 1):.0%}" if len(pats) else "-",
                          "by_patient": pats[pats > 0].to_dict()})
    ct = pd.DataFrame(cell_rows)
    log("cell level:\n" + ct.drop(columns="by_patient").to_string(index=False))
    for r in cell_rows:
        log(f"  {r['category']}: {r['by_patient']}")
    ct.to_csv(OUT / f"cells_{name}.csv", index=False)
    return t, ct


res = {}
res["stage2_harmony"] = analyse("stage2_harmony", "tirosh_malignant_annotated.h5ad", "tsoi_state", "patient",
                                ["leiden_r0.3", "leiden_r0.1", "leiden_r0.2", "leiden_r0.5", "leiden_r0.8", "leiden_r1.0"])
res["stage3_bbknn"] = analyse("stage3_bbknn", "tirosh_malignant_bbknn.h5ad", "tsoi_state_bbknn", "patient",
                              [f"leiden_bbknn_r{r}" for r in ("0.1", "0.2", "0.3", "0.5", "0.8")])
res["layer1"] = analyse("layer1_tirosh_counts", "layer1_tirosh_subset_full_annotated.h5ad", "tsoi_state", "samples",
                        ["leiden_r0.1"])
res["layer2"] = analyse("layer2_new_patients", "jerby_arnon_layer2_full_annotated.h5ad", "tsoi_state", "samples",
                        ["leiden_r0.1"])
(OUT / "log.txt").write_text("\n".join(lines) + "\n")
