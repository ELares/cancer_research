#!/usr/bin/env python3
"""Interactive corpus + simulation dashboard (#354).

A Streamlit front-end that turns the repo from a static archive into a usable
research tool: explore the corpus (filters, mechanism/cancer/evidence views, the
mechanism x cancer matrix) and run a single-cell ferroptosis parameter sweep.

Run:
    pip install -r requirements-dashboard.txt
    streamlit run scripts/dashboard.py

All aggregation logic lives in `scripts/dashboard_data.py` (stdlib-only, unit-
tested in CI). Streamlit + pandas are UI-only, optional dependencies (NOT in
requirements-lock.txt). The simulation sweep needs the compiled `ferroptosis_core`
extension; if it is not importable the tab degrades to the committed
prior-predictive intervals (read-only), per the issue's "committed outputs first,
live runs optional".
"""

import json
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dashboard_data as dd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]


@st.cache_data
def _records():
    return dd.load_index()


def corpus_tab(records):
    st.subheader("Corpus exploration")
    mech_opts = list(dd.value_counts(records, "mechanisms"))
    canc_opts = list(dd.value_counts(records, "cancer_types"))
    ev_opts = list(dd.value_counts(records, "evidence_level"))
    yrs = dd.year_histogram(records)
    ymin, ymax = (min(yrs), max(yrs)) if yrs else (2001, 2026)

    with st.sidebar:
        st.markdown("### Filters")
        f_mech = st.multiselect("Mechanism", mech_opts)
        f_canc = st.multiselect("Cancer type", canc_opts)
        f_ev = st.multiselect("Evidence level", ev_opts)
        f_year = st.slider("Year range", ymin, ymax, (ymin, ymax))

    filt = dd.filter_records(records, f_mech or None, f_canc or None, f_ev or None, f_year)
    s = dd.summary_stats(filt)
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Articles", f"{s['n_records']:,}")
    c2.metric("Mechanisms", s["n_mechanisms"])
    c3.metric("Cancer types", s["n_cancer_types"])
    c4.metric("Evidence-tagged", f"{s['n_evidence_tagged']:,}")

    st.markdown("**Mechanisms** (filtered)")
    st.bar_chart(pd.Series(dd.value_counts(filt, "mechanisms")))
    col_a, col_b = st.columns(2)
    with col_a:
        st.markdown("**Evidence tiers**")
        st.bar_chart(pd.Series(dd.value_counts(filt, "evidence_level")))
    with col_b:
        st.markdown("**Publications per year**")
        st.line_chart(pd.Series(dd.year_histogram(filt)))

    st.markdown("**Mechanism x cancer matrix** (top 10 x top 10, filtered)")
    matrix = dd.mechanism_cancer_matrix(filt, top_mech=10, top_cancer=10)
    if matrix:
        mechs = sorted({m for m, _ in matrix}, key=lambda m: -dd.value_counts(filt, "mechanisms").get(m, 0))
        cancers = sorted({c for _, c in matrix}, key=lambda c: -dd.value_counts(filt, "cancer_types").get(c, 0))
        df = pd.DataFrame(0, index=mechs, columns=cancers)
        for (m, c), n in matrix.items():
            df.loc[m, c] = n
        st.dataframe(df.style.background_gradient(cmap="Blues"), width="stretch")

    st.markdown(f"**Articles** ({len(filt):,})")
    cols = ["pmid", "year", "title", "journal", "mechanisms", "cancer_types", "evidence_level", "cited_by_count"]
    table = pd.DataFrame([{k: r.get(k) for k in cols} for r in filt])
    for lc in ("mechanisms", "cancer_types"):
        if lc in table:
            table[lc] = table[lc].apply(lambda v: ", ".join(v) if isinstance(v, list) else v)
    st.dataframe(table, width="stretch", height=400)


def _load_json(rel):
    p = REPO_ROOT / rel
    return json.loads(p.read_text()) if p.exists() else None


def simulation_tab():
    st.subheader("Single-cell ferroptosis parameter sweep")
    try:
        import ferroptosis_core as fc
        have_fc = True
    except ImportError:
        have_fc = False

    if have_fc:
        st.caption("Live `ferroptosis_core.sim_batch` sweep.")
        phenos = ["Glycolytic", "OXPHOS", "Persister"]
        treatments = ["RSL3", "SDT", "PDT", "Control"]
        pheno = st.selectbox("Phenotype", phenos)
        treat = st.selectbox("Treatment", treatments, index=0)
        ranges = _load_json("analysis/prcc-results.json")
        pr = None
        if ranges:
            for v in [ranges] + list(ranges.values() if isinstance(ranges, dict) else []):
                if isinstance(v, dict) and "parameter_ranges" in v:
                    pr = v["parameter_ranges"]
                    break
        param = st.selectbox("Swept parameter", sorted(pr) if pr else ["lp_propagation"])
        lo, hi = (pr[param] if pr and param in pr else [0.0, 1.0])
        n_pts = st.slider("Sweep points", 5, 25, 11)
        n_cells = st.select_slider("Cells per point", [1000, 2000, 4000, 8000], value=2000)
        xs = [lo + (hi - lo) * i / (n_pts - 1) for i in range(n_pts)]
        ys = [fc.sim_batch(pheno, treat, n=n_cells, seed=42, **{param: x})["death_rate"] for x in xs]
        st.line_chart(pd.DataFrame({"death_rate": ys}, index=[round(x, 4) for x in xs]))
        st.caption(f"{pheno} x {treat}: death rate vs {param} over [{lo}, {hi}] (seed 42, n={n_cells}).")
    else:
        st.info(
            "The compiled `ferroptosis_core` extension is not installed, so showing the "
            "committed prior-predictive death-rate intervals (read-only). Build the extension "
            "(see simulations/ferroptosis-python/) for the live sweep."
        )
        intervals = _load_json("analysis/uncertainty-intervals.json")
        if intervals:
            st.json(intervals)
        else:
            st.write("Committed analysis outputs are under `analysis/`; the prior-predictive "
                     "intervals are documented in `analysis/uncertainty-intervals-report.md`.")


def main():
    st.set_page_config(page_title="Cancer-research dashboard", layout="wide")
    st.title("Cancer-research corpus + simulation dashboard")
    st.caption("Issue #354. Corpus index: `corpus/INDEX.jsonl`. Read the MODEL_CARD for simulation scope/caveats.")
    records = _records()
    tab1, tab2 = st.tabs(["Corpus", "Simulation sweep"])
    with tab1:
        corpus_tab(records)
    with tab2:
        simulation_tab()


if __name__ == "__main__":
    main()
