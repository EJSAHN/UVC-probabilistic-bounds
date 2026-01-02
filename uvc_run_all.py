#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
uvc_run_all.py
Reproducible UV-C analysis + physics pipeline (GitHub-safe: no hardcoded paths)

What it does:
- Loads Coffee.xlsx and Cacao.xlsx
- Computes survival_ratio = area / mean(control area per isolate) (with epsilon)
- Defines alive_flag = survival_ratio >= tau
- Writes summary tables, ANOVA (with effect size), Coffee (~10 min) CP95 table
- Performs TAU sensitivity analysis (tau_min..tau_max, step)
- Runs leakage-safe ML (GroupKFold by plate) + dumps hyperparameters/seed
- Translates Coffee CP95 upper bound into physics bounds (sigma_LB, t90/t99/t99.9) and photonic budget

This script is designed for:
- Local reproduction (paper figures/tables can be regenerated from these outputs)
- GitHub public repository (raw Excel inputs can remain outside the repo and be referenced by path)
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, Tuple, Optional, List

import numpy as np
import pandas as pd

# stats
from statsmodels.stats.proportion import proportion_confint
import statsmodels.api as sm
import statsmodels.formula.api as smf

# ML
from sklearn.model_selection import GroupKFold
from sklearn.metrics import accuracy_score, f1_score, r2_score, mean_squared_error
from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor


ALL_FEATS = ["area", "perimeter", "length", "width", "lwr", "circularity", "is_cg"]
SHAPE_FEATS = ["lwr", "circularity", "is_cg"]  # size-free




def dedupe_columns(df: pd.DataFrame, numeric_cols: Iterable[str]) -> pd.DataFrame:
    """
    Coalesce duplicate column names created by header canonicalization.

    If a numeric column appears multiple times (duplicate names), we:
      - convert each to numeric (coerce)
      - take first non-null per row (left-to-right)
      - drop all duplicates and keep one merged column

    If a non-numeric column appears multiple times, we keep the first and drop the rest.
    """
    numeric_cols = set(numeric_cols)
    cols = list(df.columns)

    # iterate unique names preserving order
    seen_names = []
    for c in cols:
        if c not in seen_names:
            seen_names.append(c)

    for name in seen_names:
        idxs = [i for i, c in enumerate(cols) if c == name]
        if len(idxs) <= 1:
            continue

        if name in numeric_cols:
            # select all duplicate columns by boolean mask
            block = df.loc[:, [c == name for c in df.columns]]
            num = block.apply(lambda s: pd.to_numeric(s, errors="coerce"))
            merged = num.bfill(axis=1).iloc[:, 0]
            # drop all occurrences
            keep_mask = [c != name for c in df.columns]
            df = df.loc[:, keep_mask].copy()
            df[name] = merged
            cols = list(df.columns)
        else:
            # keep first occurrence only
            keep = []
            kept_one = False
            for c in df.columns:
                if c == name:
                    if not kept_one:
                        keep.append(True)
                        kept_one = True
                    else:
                        keep.append(False)
                else:
                    keep.append(True)
            df = df.loc[:, keep].copy()
            cols = list(df.columns)

    return df
# ---------------------------
# Utilities
# ---------------------------

def log(msg: str) -> None:
    print(msg, flush=True)

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p

def format_pvalue(p: float) -> str:
    """Readable p-value formatting for tables."""
    if p is None or (isinstance(p, float) and (np.isnan(p) or np.isinf(p))):
        return "NA"
    if p < 1e-300:
        return "<1e-300"
    if p < 1e-4:
        return f"{p:.2e}"
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"

def canonicalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Map likely headers to internal names."""
    colmap = {}
    for c in df.columns:
        cl = str(c).strip().lower()
        cl = re.sub(r"\s+", " ", cl)

        # strip units in [] or ()
        cl_stripped = re.sub(r"\(.*?\)|\[.*?\]", "", cl).strip()

        if "isolate" in cl_stripped:
            colmap[c] = "isolate"
        elif "plate" in cl_stripped:
            colmap[c] = "plate"
        elif "treatment" in cl_stripped:
            colmap[c] = "treatment"
        elif "area" in cl_stripped:
            colmap[c] = "area"
        elif "perimeter" in cl_stripped:
            colmap[c] = "perimeter"
        elif re.search(r"\blength\b", cl_stripped):
            colmap[c] = "length"
        elif "width" in cl_stripped:
            colmap[c] = "width"
        elif "lwr" in cl:
            colmap[c] = "lwr"
        elif "circularity" in cl:
            colmap[c] = "circularity"
        elif ("distance" in cl_stripped and "cg" in cl) or "is&cg" in cl or "is cg" in cl.replace("&", " "):
            colmap[c] = "is_cg"
        else:
            colmap[c] = c

    return df.rename(columns=colmap)

def parse_band_seconds(text: str) -> Tuple[str, Optional[float]]:
    """
    Return (band, seconds) inferred from treatment string.
    - band in {"UVC", "UVB", "Sonication", "Control", "Other"}
    - seconds numeric if present; else None
    """
    if text is None or str(text).strip() == "":
        return ("Other", None)

    t = str(text).lower()
    band = "Other"
    if "uv-c" in t or "uvc" in t:
        band = "UVC"
    elif "uv-b" in t or "uvb" in t:
        band = "UVB"
    elif "sonic" in t:
        band = "Sonication"
    elif "control" in t or "mock" in t or "untreated" in t:
        band = "Control"

    sec = None
    m = re.search(r"(\d+(?:\.\d+)?)\s*(?:m|min|mins|minute|minutes)\b", t)
    if m:
        sec = float(m.group(1)) * 60.0
    else:
        m2 = re.search(r"(\d+(?:\.\d+)?)\s*(?:s|sec|second|seconds)\b", t)
        if m2:
            sec = float(m2.group(1))

    if band == "Control":
        sec = 0.0
    return (band, sec)

def load_panel(xlsx_path: Path, pathosystem_name: str) -> pd.DataFrame:
    if not xlsx_path.exists():
        raise FileNotFoundError(f"Input not found: {xlsx_path}")

    df = pd.read_excel(xlsx_path)
    df = canonicalize_columns(df)
    # Coalesce duplicate names introduced by canonicalization (common in Excel exports)
    df = dedupe_columns(df, numeric_cols=ALL_FEATS)
    df["pathosystem"] = pathosystem_name

    band_list, sec_list = [], []
    for x in df.get("treatment", []):
        b, s = parse_band_seconds(x)
        band_list.append(b)
        sec_list.append(s)

    df["band"] = band_list
    df["seconds"] = sec_list
    df["is_control"] = df["band"].eq("Control")

    for k in ALL_FEATS:
        if k in df.columns:
            df[k] = pd.to_numeric(df[k], errors="coerce")

    # Ensure plate exists for group-safe CV
    if "plate" not in df.columns:
        df["plate"] = "NA"

    return df

def make_survival(df_all: pd.DataFrame, tau: float) -> pd.DataFrame:
    """Compute control means and survival ratio; set alive_flag using tau."""
    ctrl = (
        df_all[df_all["is_control"]]
        .groupby(["pathosystem", "isolate"], dropna=False)["area"]
        .mean()
        .rename("ctrl_area_mean")
        .reset_index()
    )
    df = df_all.merge(ctrl, on=["pathosystem", "isolate"], how="left")

    eps = 1e-9
    df["survival_ratio"] = df["area"] / (df["ctrl_area_mean"] + eps)
    df["alive_flag"] = (df["survival_ratio"] >= float(tau)).astype(int)
    df["treated_flag"] = (~df["is_control"]).astype(int)
    return df

def summary_table(df: pd.DataFrame, out_tables: Path) -> Path:
    gcols = ["pathosystem", "isolate", "treatment"]
    agg = df.groupby(gcols)[ALL_FEATS + ["survival_ratio"]].agg(["count", "mean", "std", "sem"])
    agg.columns = ["_".join([a for a in col if a]) for col in agg.columns.to_flat_index()]
    agg = agg.reset_index()
    out = out_tables / "summary_by_pathosystem_isolate_treatment.csv"
    agg.to_csv(out, index=False)
    return out

def anova_by_isolate(df: pd.DataFrame, out_tables: Path) -> Path:
    """
    One-way (control vs treated) by isolate, per trait.
    Adds effect size (eta^2) and formatted p-value column.
    """
    rows = []
    for (ps, iso), sub in df.groupby(["pathosystem", "isolate"], dropna=False):
        if sub["treatment"].nunique() < 2:
            continue

        sub = sub.copy()
        sub["grp"] = np.where(sub["is_control"], "control", "treated")

        for trait in ALL_FEATS:
            if trait not in sub.columns:
                continue
            if sub[trait].notna().sum() < 5:
                continue

            try:
                # log1p to stabilize area-like metrics
                sub["_y"] = np.log1p(sub[trait])
                fit = smf.ols("_y ~ C(grp)", data=sub).fit()
                table = sm.stats.anova_lm(fit, typ=2)

                p = float(table.loc["C(grp)", "PR(>F)"])
                ss_effect = float(table.loc["C(grp)", "sum_sq"])
                ss_resid = float(table.loc["Residual", "sum_sq"])
                eta_sq = ss_effect / (ss_effect + ss_resid) if (ss_effect + ss_resid) > 0 else np.nan

                rows.append(
                    {
                        "pathosystem": ps,
                        "isolate": iso,
                        "trait": trait,
                        "p_value": p,
                        "p_value_fmt": format_pvalue(p),
                        "eta_sq": eta_sq,
                        "n": int(sub.shape[0]),
                    }
                )
            except Exception:
                continue

    out = out_tables / "anova_by_isolate_treatment.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    return out

def coffee_cp95_from_df(df: pd.DataFrame, tau: float, uvc_seconds: float, window_s: float) -> pd.DataFrame:
    """
    Coffee, UVC only, treated, ~10 min window.
    Returns isolate-level CP95 on alive proportion using the current tau.
    """
    d0 = df[(df["pathosystem"].eq("Coffee")) & (~df["is_control"])]
    d0 = d0[d0["band"].eq("UVC")]

    # recompute alive_flag for the requested tau (df may have been built with baseline tau)
    d0 = d0.copy()
    d0["alive_flag_tau"] = (d0["survival_ratio"] >= float(tau)).astype(int)

    d1 = d0.copy()
    if d0["seconds"].notna().any():
        d1 = d0[np.abs(d0["seconds"].astype(float) - float(uvc_seconds)) <= float(window_s)]
        if d1.empty:
            d1 = d0.copy()
    # else: no seconds parsed; use all UVC treated as a fallback

    rows = []
    for iso, g in d1.groupby("isolate", dropna=False):
        n = int(g.shape[0])
        k = int(g["alive_flag_tau"].sum())
        lo, up = proportion_confint(count=k, nobs=n, alpha=0.05, method="beta")
        rows.append(
            {
                "isolate": iso,
                "tau": float(tau),
                "N": n,
                "alive": k,
                "p_hat": (k / n) if n > 0 else np.nan,
                "cp95_low": float(lo),
                "cp95_up": float(up),
            }
        )
    return pd.DataFrame(rows).sort_values(["isolate"])

def tau_sensitivity_coffee(df: pd.DataFrame, out_tables: Path, uvc_seconds: float, window_s: float,
                           tau_min: float, tau_max: float, tau_step: float, tau_baseline: float) -> Tuple[Path, Path]:
    taus = np.round(np.arange(tau_min, tau_max + 1e-12, tau_step), 6).tolist()
    long_rows = []
    for t in taus:
        tab = coffee_cp95_from_df(df, tau=t, uvc_seconds=uvc_seconds, window_s=window_s)
        long_rows.append(tab)
    long_df = pd.concat(long_rows, ignore_index=True) if long_rows else pd.DataFrame()

    out_long = out_tables / "coffee_10min_external_validation_cp95_by_tau.csv"
    long_df.to_csv(out_long, index=False)

    # robustness summary per isolate
    if not long_df.empty:
        base = long_df[long_df["tau"].eq(float(tau_baseline))][["isolate", "cp95_up"]].rename(columns={"cp95_up": "cp95_up_tau_baseline"})
        summ = (
            long_df.groupby("isolate", dropna=False)["cp95_up"]
            .agg(cp95_up_min="min", cp95_up_max="max")
            .reset_index()
        )
        summ = summ.merge(base, on="isolate", how="left")
        summ["cp95_up_range"] = summ["cp95_up_max"] - summ["cp95_up_min"]
    else:
        summ = pd.DataFrame(columns=["isolate","cp95_up_min","cp95_up_max","cp95_up_tau_baseline","cp95_up_range"])

    out_robust = out_tables / "coffee_10min_external_validation_cp95_tau_robustness.csv"
    summ.to_csv(out_robust, index=False)
    return out_long, out_robust

def ml_reports(df: pd.DataFrame, out_metrics: Path, seed: int = 42) -> Tuple[Path, Path]:
    """
    Shape-only + (optional seconds) features; GroupKFold by plate.
    Dumps metrics + model hyperparameters.
    """
    feats = SHAPE_FEATS.copy()
    dd = df.copy()
    if "seconds" in dd.columns:
        dd["seconds_filled"] = dd["seconds"].fillna(0.0)
        feats = feats + ["seconds_filled"]

    use = dd[feats + ["pathosystem", "plate", "survival_ratio"]].dropna()
    out_cls = out_metrics / "ml_classification_host_origin_GROUPSAFE.json"
    out_reg = out_metrics / "ml_regression_survival_ratio_GROUPSAFE.json"

    if use.empty:
        out_cls.write_text(json.dumps({"status": "EMPTY"}, indent=2), encoding="utf-8")
        out_reg.write_text(json.dumps({"status": "EMPTY"}, indent=2), encoding="utf-8")
        return out_cls, out_reg

    X = use[feats].copy()
    y_cls = use["pathosystem"].astype("category").cat.codes  # Coffee/Cacao
    y_reg = use["survival_ratio"].astype(float)
    groups = use["plate"].astype(str)

    n_groups = groups.nunique()
    n_splits = min(5, max(2, n_groups))
    gkf = GroupKFold(n_splits=n_splits)

    # Classifier
    accs, f1s = [], []
    clf = GradientBoostingClassifier(random_state=seed)
    clf_params = clf.get_params()
    for tr, te in gkf.split(X, y_cls, groups):
        clf = GradientBoostingClassifier(random_state=seed)
        clf.fit(X.iloc[tr], y_cls.iloc[tr])
        yp = clf.predict(X.iloc[te])
        accs.append(accuracy_score(y_cls.iloc[te], yp))
        f1s.append(f1_score(y_cls.iloc[te], yp, average="macro"))

    cls_report = {
        "features": feats,
        "cv": f"{n_splits}-fold GroupKFold by plate",
        "seed": seed,
        "model": "GradientBoostingClassifier",
        "hyperparameters": clf_params,
        "Accuracy_mean": float(np.mean(accs)),
        "Accuracy_sd": float(np.std(accs, ddof=1)) if len(accs) > 1 else 0.0,
        "MacroF1_mean": float(np.mean(f1s)),
        "MacroF1_sd": float(np.std(f1s, ddof=1)) if len(f1s) > 1 else 0.0,
    }
    out_cls.write_text(json.dumps(cls_report, indent=2), encoding="utf-8")

    # Regressor
    r2s, rmses = [], []
    reg0 = GradientBoostingRegressor(random_state=seed)
    reg_params = reg0.get_params()
    for tr, te in gkf.split(X, y_reg, groups):
        reg = GradientBoostingRegressor(random_state=seed)
        reg.fit(X.iloc[tr], y_reg.iloc[tr])
        yp = reg.predict(X.iloc[te])
        r2s.append(r2_score(y_reg.iloc[te], yp))
        rmses.append(math.sqrt(mean_squared_error(y_reg.iloc[te], yp)))

    reg_report = {
        "features": feats,
        "cv": f"{n_splits}-fold GroupKFold by plate",
        "seed": seed,
        "model": "GradientBoostingRegressor",
        "hyperparameters": reg_params,
        "R2_mean": float(np.mean(r2s)),
        "R2_sd": float(np.std(r2s, ddof=1)) if len(r2s) > 1 else 0.0,
        "RMSE_mean": float(np.mean(rmses)),
        "RMSE_sd": float(np.std(rmses, ddof=1)) if len(rmses) > 1 else 0.0,
    }
    out_reg.write_text(json.dumps(reg_report, indent=2), encoding="utf-8")

    return out_cls, out_reg

def physics_from_cp95(cp95_df: pd.DataFrame,
                      out_tables: Path,
                      out_metrics: Path,
                      intensity_mw_cm2: float,
                      wavelength_nm: float,
                      uvc_seconds: float,
                      tau: float) -> Tuple[Path, Path]:
    """
    Physics translation based on CP95 upper bound on survival probability.
    Model: S = exp(-sigma * H), H = I * t
    We compute sigma_LB = -ln(S_up) / H, then t90/t99/t99.9 at the same irradiance.
    """
    I = float(intensity_mw_cm2) * 1e-3  # W/cm^2
    t0 = float(uvc_seconds)
    H = I * t0  # J/cm^2

    # photon budget
    h = 6.62607015e-34
    c = 2.99792458e8
    lam = float(wavelength_nm) * 1e-9
    E_photon = (h * c) / lam  # J/photon
    N_photons = H / E_photon  # photons/cm^2

    rows = []
    for _, r in cp95_df.iterrows():
        iso = str(r["isolate"])
        s_up = float(r["cp95_up"])
        s_up = min(max(s_up, 1e-12), 1.0)

        sigma_lb = 0.0 if s_up >= 1.0 else (-math.log(s_up) / H)  # cm^2/J

        def t_for_kill(frac: float) -> float:
            s_target = 1.0 - frac
            if sigma_lb <= 0:
                return float("inf")
            return (-math.log(s_target)) / (sigma_lb * I)

        rows.append({
            "isolate": iso,
            "tau": float(tau),
            "cp95_up": s_up,
            "sigma_LB_cm2_per_J": sigma_lb,
            "t90_s": t_for_kill(0.90),
            "t99_s": t_for_kill(0.99),
            "t999_s": t_for_kill(0.999),
        })

    tab = pd.DataFrame(rows).sort_values("isolate")
    out_tab = out_tables / "physics_bounds_by_isolate.csv"
    tab.to_csv(out_tab, index=False)

    summary = {
        "intensity_mw_cm2": float(intensity_mw_cm2),
        "wavelength_nm": float(wavelength_nm),
        "uvc_seconds": float(uvc_seconds),
        "H_J_cm2": float(H),
        "photon_fluence_photons_per_cm2": float(N_photons),
        "tau_baseline": float(tau),
        "note": "If cp95_up==1.0 then sigma_LB=0 and times are infinite (no bounded guarantee under this model)."
    }
    out_summary = out_metrics / "physics_summary.json"
    out_summary.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return out_tab, out_summary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=str, default="", help="Folder containing Coffee.xlsx and Cacao.xlsx (recommended allowing raw data to live outside the repo).")
    ap.add_argument("--coffee", type=str, default="", help="Path to Coffee.xlsx (overrides --data-dir).")
    ap.add_argument("--cacao", type=str, default="", help="Path to Cacao.xlsx (overrides --data-dir).")
    ap.add_argument("--out-dir", type=str, default="outputs", help="Output directory.")
    ap.add_argument("--tau", type=float, default=0.05, help="Alive threshold: survival_ratio >= tau.")
    ap.add_argument("--tau-min", type=float, default=0.01, help="Tau sensitivity minimum.")
    ap.add_argument("--tau-max", type=float, default=0.10, help="Tau sensitivity maximum.")
    ap.add_argument("--tau-step", type=float, default=0.01, help="Tau sensitivity step.")
    ap.add_argument("--uvc-seconds", type=float, default=600.0, help="Coffee ~10 min window center in seconds.")
    ap.add_argument("--uvc-window-seconds", type=float, default=60.0, help="Half-window size in seconds (±).")
    ap.add_argument("--uvc-intensity-mw-cm2", type=float, default=0.58, help="Irradiance in mW/cm^2.")
    ap.add_argument("--uvc-wavelength-nm", type=float, default=275.0, help="Wavelength in nm.")
    ap.add_argument("--seed", type=int, default=42, help="Random seed for ML.")
    args = ap.parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_tables = ensure_dir(out_dir / "tables")
    out_metrics = ensure_dir(out_dir / "metrics")
    out_logs = ensure_dir(out_dir / "logs")

    # Resolve inputs
    coffee_path = Path(args.coffee) if args.coffee else None
    cacao_path = Path(args.cacao) if args.cacao else None
    if args.data_dir:
        dd = Path(args.data_dir)
        if coffee_path is None or not str(coffee_path):
            coffee_path = dd / "Coffee.xlsx"
        if cacao_path is None or not str(cacao_path):
            cacao_path = dd / "Cacao.xlsx"

    if not coffee_path or not cacao_path:
        log("[ERROR] Provide inputs via --data-dir or --coffee/--cacao.")
        return 2

    coffee_path = coffee_path.resolve()
    cacao_path = cacao_path.resolve()

    log(f"[INFO] Coffee = {coffee_path}")
    log(f"[INFO] Cacao  = {cacao_path}")
    log(f"[INFO] Out    = {out_dir}")
    log(f"[INFO] tau    = {args.tau}")

    # Load
    cacao = load_panel(cacao_path, "Cacao")
    coffee = load_panel(coffee_path, "Coffee")
    df0 = pd.concat([cacao, coffee], ignore_index=True)

    # Compute survival at baseline tau
    df = make_survival(df0, tau=float(args.tau))

    # Tables
    summary_path = summary_table(df, out_tables)
    log(f"[OK] {summary_path.name}")

    anova_path = anova_by_isolate(df, out_tables)
    log(f"[OK] {anova_path.name}")

    # Coffee CP95 baseline tau
    cp = coffee_cp95_from_df(df, tau=float(args.tau), uvc_seconds=float(args.uvc_seconds), window_s=float(args.uvc_window_seconds))
    out_cp = out_tables / "coffee_10min_external_validation_cp95.csv"
    cp.to_csv(out_cp, index=False)
    log(f"[OK] {out_cp.name}")

    # Tau sensitivity
    out_long, out_rob = tau_sensitivity_coffee(
        df,
        out_tables=out_tables,
        uvc_seconds=float(args.uvc_seconds),
        window_s=float(args.uvc_window_seconds),
        tau_min=float(args.tau_min),
        tau_max=float(args.tau_max),
        tau_step=float(args.tau_step),
        tau_baseline=float(args.tau),
    )
    log(f"[OK] {out_long.name}")
    log(f"[OK] {out_rob.name}")

    # ML reports
    out_cls, out_reg = ml_reports(df, out_metrics, seed=int(args.seed))
    log(f"[OK] {out_cls.name}")
    log(f"[OK] {out_reg.name}")

    # Physics translation from baseline tau CP95
    out_phys_tab, out_phys_summary = physics_from_cp95(
        cp95_df=cp,
        out_tables=out_tables,
        out_metrics=out_metrics,
        intensity_mw_cm2=float(args.uvc_intensity_mw_cm2),
        wavelength_nm=float(args.uvc_wavelength_nm),
        uvc_seconds=float(args.uvc_seconds),
        tau=float(args.tau),
    )
    log(f"[OK] {out_phys_tab.name}")
    log(f"[OK] {out_phys_summary.name}")

    # Full row-level export
    out_full = out_tables / "summary_full_dataset_with_survival.csv"
    df.to_csv(out_full, index=False)
    log(f"[OK] {out_full.name}")

    # Run metadata
    run_info = {
        "timestamp_utc": datetime.utcnow().isoformat() + "Z",
        "inputs": {"coffee": str(coffee_path), "cacao": str(cacao_path)},
        "outputs": {"out_dir": str(out_dir)},
        "params": {
            "tau": float(args.tau),
            "tau_min": float(args.tau_min),
            "tau_max": float(args.tau_max),
            "tau_step": float(args.tau_step),
            "uvc_seconds": float(args.uvc_seconds),
            "uvc_window_seconds": float(args.uvc_window_seconds),
            "uvc_intensity_mw_cm2": float(args.uvc_intensity_mw_cm2),
            "uvc_wavelength_nm": float(args.uvc_wavelength_nm),
            "seed": int(args.seed),
        },
        "notes": [
            "This pipeline is GitHub-safe: raw Excel inputs can be kept outside the repository.",
            "Reviewer-driven updates included: tau sensitivity analysis table(s), ANOVA effect size (eta^2), ML hyperparameter dump, and physics bounds.",
        ],
    }
    (out_logs / "RUN_INFO.json").write_text(json.dumps(run_info, indent=2), encoding="utf-8")
    log("[DONE] Pipeline completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
