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
        # use_container_width (not width="stretch") for cross-version compatibility:
        # the stlite/Pyodide demo (#565) bundles Streamlit 1.39, where the string
        # `width` API does not exist; use_container_width works 1.39->current
        # (deprecation-warned, not an error, in the newest pinned local Streamlit).
        # Keep numeric counts readable without pandas Styler's matplotlib
        # dependency: the browser's bundled versions disagree on get_cmap.
        st.dataframe(df, use_container_width=True)

    st.markdown(f"**Articles** ({len(filt):,})")
    cols = ["pmid", "year", "title", "journal", "mechanisms", "cancer_types", "evidence_level", "cited_by_count"]
    table = pd.DataFrame([{k: r.get(k) for k in cols} for r in filt])
    for lc in ("mechanisms", "cancer_types"):
        if lc in table:
            table[lc] = table[lc].apply(lambda v: ", ".join(v) if isinstance(v, list) else v)
    st.dataframe(table, use_container_width=True, height=400)


@st.cache_data
def _census():
    return dd.load_census()


def census_tab():
    """The census, at the only resolution a browser can hold: aggregates.

    Record-level browsing is not offered and its absence is stated rather than
    worked around. 5,187,265 records cannot be loaded client-side, and the
    census is gitignored besides -- so a tab that appeared to browse it would
    be browsing something else.
    """
    st.subheader("Census")
    c = _census()
    head = dd.census_headline(c.get("design"))
    if head is None:
        st.warning(
            "The committed census aggregates are missing from `analysis/`. "
            "Regenerate them with the `scripts/census_*.py` generators; this "
            "panel shows nothing rather than showing a partial census that "
            "looks complete."
        )
        return

    a, b, d = st.columns(3)
    a.metric("Indexed census records", f"{head['census']:,}")
    b.metric("Trial-labelled records", f"{head['trials']:,}")
    d.metric("Undetermined design", f"{head['undetermined']:,}")
    st.caption(
        f"The trial share has two denominators and both are shown because "
        f"either alone misleads in a predictable direction: "
        f"**{head['share_of_census']}%** of the whole census, "
        f"**{head['share_of_classifiable']}%** of the "
        f"{head['classifiable']:,} records carrying a design-informative "
        f"label at all. Study design is read from NLM publication types and "
        f"MeSH check tags. Counts describe labelled publications, not distinct "
        f"clinical trials or treatment effectiveness."
    )

    rows = dd.census_mechanism_rows(c.get("profile"))
    if rows:
        st.markdown("#### Mechanisms, ordered by clinical-trial share")
        st.caption(
            "Trial share is the proportion of records carrying a trial label "
            "within each mapped mechanism. Descriptor breadth and overlapping "
            "tags still affect these comparisons. Neither trial share nor "
            "publication volume ranks treatment effectiveness."
        )
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
        st.caption(
            "TTFields and bioelectric modulation are absent from this "
            "descriptor mapping. They are unmeasurable here, not zero research."
        )

    growth = c.get("growth")
    if growth and growth.get("union_growth"):
        st.markdown("#### Growth, against the denominator a growth claim needs")
        st.caption(
            f"Cancer literature as a whole grew x{growth['field_growth']} "
            f"between {growth['start_year']} and {growth['end_year']}; the "
            f"mechanisms tracked here grew x{growth['union_growth']}, which is "
            f"x{growth['mechanisms_over_field']} the field. A corpus of "
            f"emerging-therapy papers outgrows all of cancer research whether "
            f"or not anything unusual happened, so the field is the wrong "
            f"comparator and the mechanisms' own rate is the right one."
        )

    st.caption(
        "These aggregate panels summarize 4,403,994 records indexed under neoplasms "
        "and adjacent descriptors, within a broader 5,187,265-record archive. "
        "Record-level browsing of the "
        "census is not offered. The Corpus tab browses the separate 4,830-record "
        "historical archive, retained as a method-comparison arm."
    )


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
        st.caption(
            "Live `ferroptosis_core.sim_batch` sweep. These four arms are what "
            "the PYTHON BINDING accepts, not what the engine has -- `sim_cell` "
            "takes no dose, depth or alpha/beta, so the arms that need them "
            "(radiation, ablation, ADC, adoptive cell, oncolytic) are "
            "deliberately unreachable here and run from the Rust binaries."
        )
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
        st.markdown(
            "**Historical prior-predictive analysis.** This report predates the "
            "current calibration studies; its statements about unavailable "
            "calibration data describe that earlier stage. These intervals are "
            "not a calibrated posterior. See the [current calibration ledger]"
            "(https://github.com/ELares/cancer_research/blob/main/"
            "simulations/calibration/CALIBRATION_STATUS.md) for the latest "
            "fits, sampling checks and remaining validation limits."
        )
        intervals = _load_json("analysis/uncertainty-intervals.json")
        if intervals:
            st.json(intervals)
        else:
            report = REPO_ROOT / "analysis/uncertainty-intervals-report.md"
            if report.exists():
                st.markdown(report.read_text(encoding="utf-8"))
            else:
                st.warning("The committed prior-predictive report is unavailable.")


def main():
    st.set_page_config(page_title="Cancer-research dashboard", layout="wide")
    st.title("Cancer-research census + simulation dashboard")
    st.caption("Census aggregates: `analysis/census-*.json`. Corpus index: "
               "`corpus/INDEX.jsonl`. Read the MODEL_CARD for simulation "
               "scope/caveats.")
    tab0, tab1, tab2 = st.tabs(["Census", "Corpus (control arm)", "Simulation sweep"])
    with tab0:
        census_tab()
    with tab1:
        corpus_tab(_records())
    with tab2:
        simulation_tab()
    st.caption("All dashboard panels are ready. Article filters apply to the historical archive in the Corpus tab.")


if __name__ == "__main__":
    main()
