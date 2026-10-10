# ==============================================================================
# scripts/04_score_cohort.py
# Purpose: compute the Resistance Alteration Score (RAS) on your own hotspot-panel cohort,
#         and associate it with clinical outcome (PFS/OS) -> real-world external validation of the public signature
# Input: cohort.csv  (absolute path; fields below)
# Dependencies: pandas, optional lifelines (survival); install into an isolated venv, do not pollute the system
# ==============================================================================
import pandas as pd
import os
from pathlib import Path

# ---------- configuration (edit these three paths to run on your own cohort) ----------
COHORT_CSV = os.path.join(str(Path(__file__).resolve().parents[1]), "example", "cohort_example.csv")
GENE_CSV   = os.path.join(str(Path(__file__).resolve().parents[1]), "example", "resistance_genes.csv")
OUT_CSV    = os.path.join(str(Path(__file__).resolve().parents[1]), "example", "expected_output", "cohort_ras_scored.csv")

# expected cohort.csv columns:
# patient_id, gene, alteration(mutation/deletion/amplification), vaf, subtype(HR+/HER2-/HER2+/TNBC),
# line(treatment line), prior_AI_months(months of prior AI exposure), pfs_months, os_months, event_pfs, event_os
# one row = one alteration of one patient; a patient may have multiple rows

def main():
    genes = pd.read_csv(GENE_CSV)
    panel = set(genes["gene"].dropna().astype(str))
    df = pd.read_csv(COHORT_CSV)

    # keep only deleterious alterations belonging to the resistance panel
    df["in_panel"] = df["gene"].isin(panel)
    df_panel = df[df["in_panel"]].copy()

    # RAS: per-patient count of resistance alterations + whether >=1
    ras = (df_panel.groupby("patient_id")
                  .size()
                  .rename("ras_n")
                  .reset_index())
    ras["ras_pos"] = (ras["ras_n"] >= 1).map({True: "RAS+", False: "RAS-"})

    # merge clinical (one row per patient)
    clin = (df[["patient_id","subtype","line","prior_AI_months",
               "pfs_months","os_months","event_pfs","event_os"]]
            .drop_duplicates("patient_id"))
    out = clin.merge(ras, on="patient_id", how="left")
    out["ras_n"]   = out["ras_n"].fillna(0)
    out["ras_pos"] = out["ras_pos"].fillna("RAS-")

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"[ok] written {OUT_CSV}")
    print(out["ras_pos"].value_counts())
    print(out.groupby("ras_pos")[["pfs_months","os_months"]].median())

    # ---- survival analysis (if lifelines is installed) ----
    try:
        from lifelines import KaplanMeierFitter, CoxPHFitter
        km = KaplanMeierFitter()
        for grp, g in out.groupby("ras_pos"):
            km.fit(g["pfs_months"], event_observed=g["event_pfs"], label=grp)
            print(f"[km] {grp}: median PFS = {km.median_survival_time_:.1f} months")
        cph = CoxPHFitter()
        try:
            cph.fit(out, duration_col="pfs_months", event_col="event_pfs",
                    formula="ras_pos + C(subtype) + line")
        except Exception:
            # small cohorts: fall back to the unadjusted univariate model
            cph.fit(out[["ras_pos", "pfs_months", "event_pfs"]],
                    duration_col="pfs_months", event_col="event_pfs")
        print(cph.summary)
    except ImportError:
        print("[note] lifelines not installed -> grouped descriptions only; run the survival analysis in R")

if __name__ == "__main__":
    main()
