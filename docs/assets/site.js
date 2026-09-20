/* The atlas reads one local, verified snapshot; no runtime or analytics downloads. */
(() => {
  "use strict";
  const $ = (selector) => document.querySelector(selector);
  const $$ = (selector) => [...document.querySelectorAll(selector)];
  const number = (value) => new Intl.NumberFormat("en-US").format(value);
  const escape = (value) =>
    String(value).replace(
      /[&<>"']/g,
      (char) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[char],
    );
  const normalized = (value) =>
    value
      .toLowerCase()
      .normalize("NFKD")
      .replace(/[\u0300-\u036f]/g, "")
      .replace(/[–—_-]/g, " ");
  let snapshot = null;
  let selected = null;
  let expanded = false;
  let geometry = "rotated_box";
  let loading = false;
  // Restore the whole unloaded view if parsing or rendering fails partway through.
  const initialPanels = $$(
    "[data-stat], #hero-pass-count, #hero-trials, #mechanism-count, #mechanism-list, #mechanism-detail, #sampling-summary, #challenge-summary, #geometry-visual, #challenge-detail, #snapshot-date, #source-list",
  ).map((node) => [node, node.innerHTML]);
  const initialDate = $("#snapshot-date").dateTime;

  function stat(name, value) {
    $$(`[data-stat="${name}"]`).forEach((node) => {
      node.textContent = value;
    });
  }

  function renderOverview() {
    const c = snapshot.census;
    stat("collection-total", number(c.collection_total));
    stat("indexed", number(c.indexed));
    stat("recovered", number(c.recovered));
    stat("mechanisms", number(snapshot.mechanisms.length));
    stat("trial-count", number(c.trials));
    stat("trial-census-share", `${c.share_of_census.toFixed(2)}%`);
    stat("trial-classifiable-share", `${c.share_of_classifiable.toFixed(2)}%`);
    stat("classifiable", number(c.classifiable));
    const challenge = snapshot.challenge;
    $("#hero-pass-count").textContent = `${challenge.passed_runs}/${challenge.total_runs}`;
    $("#hero-trials").innerHTML = challenge.groups
      .map(
        (group) =>
          `<div class="trial-row"><span class="trial-label">${escape(group.label)}</span><span class="trial-marks" aria-hidden="true">${group.runs.map((run) => `<span class="trial-mark ${run.passed ? "pass" : "fail"}">${run.passed ? "✓" : "!"}</span>`).join("")}</span><span class="trial-score" aria-label="${group.passed_runs} of ${group.total_runs} runs passed">${group.passed_runs}/${group.total_runs}</span></div>`,
      )
      .join("");
    const runs = snapshot.sampling.runs;
    const ess = runs.map((run) => run.ess);
    $("#sampling-summary").innerHTML =
      `<strong>${runs.filter((run) => run.passed).length} / ${runs.length} runs passed</strong><span>Weight ESS ${Math.min(...ess).toFixed(1)}–${Math.max(...ess).toFixed(1)} · ${number(runs[0].attempts)} attempts per run</span>`;
    $("#challenge-summary").innerHTML =
      `<strong>${challenge.passed_runs} / ${challenge.total_runs} runs passed</strong><span>All runs had to pass the fixed overall rule</span>`;
    const date = new Date(`${snapshot.snapshot_date}T12:00:00Z`);
    $("#snapshot-date").dateTime = snapshot.snapshot_date;
    $("#snapshot-date").textContent = new Intl.DateTimeFormat("en-GB", {
      day: "numeric",
      month: "long",
      year: "numeric",
      timeZone: "UTC",
    }).format(date);
    $("#source-list").innerHTML = snapshot.sources
      .map(
        (source) =>
          `<li><a href="${escape(source.url)}">${escape(source.path)} ↗</a><code>SHA256 ${escape(source.sha256)}</code></li>`,
      )
      .join("");
  }

  function filteredMechanisms() {
    const query = normalized($("#mechanism-search").value.trim());
    const sort = $("#mechanism-sort").value;
    const rows = snapshot.mechanisms.filter((row) => normalized(`${row.label} ${row.id}`).includes(query));
    rows.sort((a, b) =>
      sort === "name"
        ? a.label.localeCompare(b.label, "en")
        : b[sort] - a[sort] || a.label.localeCompare(b.label, "en"),
    );
    return rows;
  }

  function renderMechanisms() {
    if (!snapshot) return;
    const rows = filteredMechanisms();
    const shown = expanded ? rows : rows.slice(0, 8);
    $("#mechanism-count").textContent = rows.length
      ? `Showing ${shown.length} of ${rows.length} mechanisms`
      : "No matching mechanisms";
    $("#show-all").hidden = rows.length <= 8;
    $("#show-all").innerHTML = expanded
      ? 'Show fewer mechanisms <span aria-hidden="true">↑</span>'
      : `Show all ${rows.length} mechanisms <span aria-hidden="true">↓</span>`;
    $("#sort-note").textContent =
      $("#mechanism-sort").value === "articles"
        ? "Article counts are not a ranking of importance or efficacy: descriptor widths differ. Trial-share bars use a common scale."
        : "Trial share = trial-labelled records / all records tagged with that mechanism. Tags can overlap.";
    if (!rows.length) {
      $("#mechanism-list").innerHTML =
        '<div class="empty-state"><h3>No mechanisms found.</h3><p>Try a broader term, such as “immunotherapy” or “antibody”.</p><button type="button" id="clear-search">Clear search</button></div>';
      $("#mechanism-detail").innerHTML =
        '<div class="panel-placeholder"><span class="eyebrow">KEEP EXPLORING</span><h3>A missing match is not a research gap.</h3><p>This atlas has a limited descriptor mapping. Try another term or inspect the source report.</p><a class="text-link" href="https://github.com/ELares/cancer_research/blob/main/analysis/census-mechanism-profile.md">Read the methods ↗</a></div>';
      $("#clear-search").addEventListener("click", () => {
        $("#mechanism-search").value = "";
        expanded = false;
        renderMechanisms();
        $("#mechanism-search").focus();
      });
      selected = null;
      return;
    }
    if (!shown.some((row) => row.id === selected)) selected = shown[0].id;
    const ceiling = Math.max(10, Math.ceil(Math.max(...snapshot.mechanisms.map((row) => row.trial_share)) / 5) * 5);
    $(".list-caption>span:last-child").textContent = `TRIAL SHARE · 0–${ceiling}%`;
    $("#mechanism-list").innerHTML = shown
      .map(
        (row, index) =>
          `<button type="button" class="mechanism-row" data-mechanism="${escape(row.id)}" aria-pressed="${row.id === selected}" aria-controls="mechanism-detail"><span class="mechanism-index" aria-hidden="true">${String(index + 1).padStart(2, "0")}</span><span><span class="mechanism-name">${escape(row.label)}</span><span class="mechanism-count">${number(row.articles)} records · ${number(row.trials)} trial-labelled</span></span><span class="share-column">${row.trial_share.toFixed(2)}%<span class="share-bar" aria-hidden="true"><i style="width:${(row.trial_share / ceiling) * 100}%"></i></span></span></button>`,
      )
      .join("");
    $$(".mechanism-row").forEach((button) =>
      button.addEventListener("click", () => {
        selected = button.dataset.mechanism;
        $$(".mechanism-row").forEach((row) =>
          row.setAttribute("aria-pressed", String(row.dataset.mechanism === selected)),
        );
        renderMechanismDetail();
        $("#selection-announcement").textContent =
          `Showing details for ${snapshot.mechanisms.find((row) => row.id === selected).label}.`;
        const heading = $("#detail-heading");
        const headingBounds = heading.getBoundingClientRect();
        const headerBottom = $(".site-header").getBoundingClientRect().bottom;
        if (
          window.matchMedia("(max-width: 600px)").matches ||
          headingBounds.top < headerBottom ||
          headingBounds.bottom > window.innerHeight
        ) {
          heading.focus({ preventScroll: true });
          $("#mechanism-detail").scrollIntoView({
            behavior: window.matchMedia("(prefers-reduced-motion: reduce)").matches ? "instant" : "smooth",
            block: "start",
          });
        }
      }),
    );
    renderMechanismDetail();
  }

  function renderMechanismDetail() {
    const row = snapshot.mechanisms.find((item) => item.id === selected);
    const period = snapshot.growth_period;
    const lookup = new Map(row.years);
    const years = Array.from({ length: period.end_year - period.start_year + 1 }, (_, i) => [
      period.start_year + i,
      lookup.get(period.start_year + i) || 0,
    ]);
    const max = Math.max(1, ...years.map((item) => item[1]));
    const growth = row.growth === null ? "Insufficient growth baseline" : `${row.growth.toFixed(2)}× publication count`;
    const sites = row.top_sites.slice(0, 3);
    $("#mechanism-detail").innerHTML = `
      <div class="detail-topline"><span class="eyebrow">MECHANISM IN FOCUS</span><span>Indexed census</span></div>
      <h3 id="detail-heading" class="detail-title" tabindex="-1">${escape(row.label)}</h3>
      <div class="detail-numbers"><div><strong>${row.trial_share.toFixed(2)}%</strong><span>clinical-trial share</span></div><div><strong>${number(row.articles)}</strong><span>tagged records · ${number(row.trials)} trial-labelled records</span></div></div>
      ${row.note ? `<p class="descriptor-note">${escape(row.note)}</p>` : ""}
      <div class="detail-section-heading"><span>Publication trend</span><span>${period.start_year}–${period.end_year}</span></div>
      <div class="trend-chart" role="img" aria-label="Annual publication counts from ${period.start_year} to ${period.end_year}. Exact values are in the table below.">${years.map(([year, count]) => `<span class="trend-column" title="${year}: ${number(count)} records"><span class="trend-bar" style="height:${(count / max) * 100}%"></span></span>`).join("")}</div>
      <div class="trend-labels"><span>${period.start_year}</span><span>${escape(growth)}</span><span>${period.end_year}</span></div>
      <details class="annual-data"><summary>View annual counts</summary><table><caption class="sr-only">${escape(row.label)} publication counts</caption><thead><tr><th scope="col">Year</th><th scope="col">Records</th></tr></thead><tbody>${years.map(([year, count]) => `<tr><th scope="row">${year}</th><td>${number(count)}</td></tr>`).join("")}</tbody></table></details>
      <div class="detail-section-heading"><span>Strongest cancer-site enrichment</span><span>vs. census site baseline</span></div>
      <ul class="site-list">${sites.length ? sites.map((site) => `<li><span>${escape(site.label)}</span><span>${site.enrichment.toFixed(2)}× · ${number(site.articles)} records</span></li>`).join("") : "<li>No site-enrichment data in this snapshot.</li>"}</ul>
      <p class="detail-footnote">Counts describe publications, not patients. Site tags overlap; enrichment is relative representation, not treatment benefit. <a href="https://github.com/ELares/cancer_research/blob/main/analysis/census-mechanism-profile.md">Check the source ↗</a></p>`;
  }

  const geometryInfo = {
    rotated_box: {
      shape:
        '<g transform="rotate(45 100 65)"><rect x="36" y="59" width="128" height="12" fill="#85ad95" fill-opacity=".38" stroke="#006b5b" stroke-width="1.5"/><path d="M36 65h128M100 59v12" stroke="#006b5b" stroke-width=".8"/></g>',
      description:
        "A narrow, rotated target with strong coordinate correlation. Two runs fail the maximum-weight limit; one also fails regional-mass accuracy. All four regions are observed.",
    },
    annular_cylinder: {
      shape:
        '<circle cx="100" cy="65" r="37.5" fill="none" stroke="#85ad95" stroke-opacity=".4" stroke-width="15"/><circle cx="100" cy="65" r="45" fill="none" stroke="#006b5b" stroke-width="1.5"/><circle cx="100" cy="65" r="30" fill="none" stroke="#006b5b" stroke-width="1.5"/>',
      description:
        "A curved annular target with eight angular regions. All three runs pass the fixed checks. This does not establish the entire within-region distribution.",
    },
    unequal_balls: {
      shape:
        '<circle cx="59" cy="82" r="25" fill="#85ad95" fill-opacity=".38" stroke="#006b5b" stroke-width="1.5"/><circle cx="143" cy="48" r="31.25" fill="#85ad95" fill-opacity=".38" stroke="#006b5b" stroke-width="1.5"/>',
      description:
        "Two disconnected seven-dimensional balls with unequal masses. All three runs pass. Finite moment checks still cannot certify the full angular distribution.",
    },
  };

  function renderChallenge() {
    if (!snapshot) return;
    const group = snapshot.challenge.groups.find((item) => item.id === geometry);
    const info = geometryInfo[geometry];
    const gates = snapshot.challenge.thresholds;
    $$(".geometry-buttons button").forEach((button) =>
      button.setAttribute("aria-pressed", String(button.dataset.geometry === geometry)),
    );
    $("#geometry-visual").innerHTML =
      `<svg viewBox="0 0 200 130" xmlns="http://www.w3.org/2000/svg"><path d="M15 15h170v100H15Z M15 65h170 M100 15v100" fill="none" stroke="#c9d5c4" stroke-width=".8"/>${info.shape}</svg><p class="schematic-label">2D schematic of a 7D target<br>Shape only; not sampled particles</p>`;
    $("#challenge-detail").innerHTML = `
      <div class="challenge-detail-heading"><h4>${escape(group.label)}</h4><span class="tag ${group.passed_runs === group.total_runs ? "tag-positive" : "tag-caution"}">${group.passed_runs}/${group.total_runs} passed</span></div>
      <div class="table-scroll" role="region" aria-label="${escape(group.label)} run results" tabindex="0"><table class="challenge-table"><caption class="sr-only">Independent ${escape(group.label)} runs. ESS means effective sample size of importance weights.</caption><thead><tr><th scope="col">Seed</th><th scope="col"><abbr title="Effective sample size of importance weights">ESS</abbr></th><th scope="col">Max weight</th><th scope="col">Result</th></tr></thead><tbody>${group.runs.map((run) => `<tr><th scope="row">${run.seed}</th><td>${run.ess.toFixed(1)}</td><td data-failed="${run.max_weight > gates.maximum_normalized_weight}">${(run.max_weight * 100).toFixed(2)}%</td><td>${run.passed ? "Pass" : "Fail"}</td></tr>`).join("")}</tbody></table></div>
      <p class="challenge-legend">Fixed limits: ESS ≥ ${gates.minimum_ess}; max weight ≤ ${(gates.maximum_normalized_weight * 100).toFixed(0)}%. <a href="https://github.com/ELares/cancer_research/blob/main/docs/COVERAGE_CHALLENGE_PLAN.md">All criteria ↗</a></p>
      <p class="challenge-description">${escape(info.description)}</p>
      <details class="annual-data"><summary>Exact regional errors & control results</summary><ul class="site-list">${group.runs.map((run) => `<li><span>Seed ${run.seed}</span><span>${run.region_error.toFixed(4)} regional error</span></li>`).join("")}</ul><p class="challenge-description">Maximum regional-mass error ≤ ${gates.maximum_region_mass_error}. All ${snapshot.challenge.controls.oracle_total} full-target controls pass. All ${snapshot.challenge.controls.negative_total} controls with deliberately missing support pass weight screens but fail truth checks.</p></details>`;
  }

  async function loadSnapshot() {
    if (loading) return;
    loading = true;
    $("#retry-data").disabled = true;
    const controller = new AbortController();
    const timeout = setTimeout(() => controller.abort(), 10000);
    try {
      const response = await fetch(new URL("assets/research-data.json", document.baseURI), {
        signal: controller.signal,
      });
      if (!response.ok) throw new Error("Snapshot unavailable");
      const data = await response.json();
      if (
        data.schema_version !== 1 ||
        !Array.isArray(data.mechanisms) ||
        !data.mechanisms.length ||
        !data.challenge?.groups?.length ||
        !data.sources?.length
      )
        throw new Error("Snapshot format invalid");
      snapshot = data;
      renderOverview();
      renderMechanisms();
      renderChallenge();
      $("#mechanism-search").disabled = false;
      $("#mechanism-sort").disabled = false;
      $$(".geometry-buttons button").forEach((button) => {
        button.disabled = false;
      });
      $("#data-error").hidden = true;
    } catch (_) {
      snapshot = null;
      selected = null;
      initialPanels.forEach(([node, markup]) => {
        node.innerHTML = markup;
      });
      $("#snapshot-date").dateTime = initialDate;
      $("#mechanism-search").disabled = true;
      $("#mechanism-sort").disabled = true;
      $$(".geometry-buttons button").forEach((button) => {
        button.disabled = true;
      });
      $("#show-all").hidden = true;
      $("#data-error").hidden = false;
      $("#hero-trials").innerHTML =
        '<p class="loading-copy">Snapshot unavailable. The source reports remain available.</p>';
      $("#mechanism-count").textContent = "Snapshot unavailable";
      $("#mechanism-detail").innerHTML =
        '<div class="panel-placeholder"><h3>The interactive snapshot is unavailable.</h3><p>Use Try again above, or follow the methods and source reports for the complete tables.</p></div>';
      $("#challenge-detail").innerHTML =
        '<p class="challenge-description">The interactive snapshot is unavailable. <a href="https://github.com/ELares/cancer_research/blob/main/analysis/calibration/proposal-coverage-challenges.md">Read all results in the source report ↗</a></p>';
      $("#sampling-summary").textContent = "See the linked report for the recorded results.";
      $("#challenge-summary").textContent = "See the linked report for the recorded results.";
    } finally {
      clearTimeout(timeout);
      loading = false;
      $("#retry-data").disabled = false;
    }
  }

  $("#mechanism-search").addEventListener("input", () => {
    expanded = false;
    renderMechanisms();
  });
  $("#mechanism-sort").addEventListener("change", renderMechanisms);
  $("#show-all").addEventListener("click", () => {
    expanded = !expanded;
    renderMechanisms();
  });
  $$(".geometry-buttons button").forEach((button) =>
    button.addEventListener("click", () => {
      geometry = button.dataset.geometry;
      renderChallenge();
    }),
  );
  $("#retry-data").addEventListener("click", loadSnapshot);
  const dialog = $("#sources-dialog");
  $("#open-sources").addEventListener("click", () => dialog.showModal());
  $("#close-sources").addEventListener("click", () => dialog.close());
  dialog.addEventListener("click", (event) => {
    if (event.target === dialog) {
      const rect = dialog.getBoundingClientRect();
      if (
        event.clientX < rect.left ||
        event.clientX > rect.right ||
        event.clientY < rect.top ||
        event.clientY > rect.bottom
      )
        dialog.close();
    }
  });
  const observer = new IntersectionObserver(
    (entries) => {
      entries.forEach((entry) => {
        if (entry.isIntersecting) {
          $$(".main-nav a").forEach((link) => {
            if (link.hash === `#${entry.target.id}`) link.setAttribute("aria-current", "location");
            else link.removeAttribute("aria-current");
          });
        }
      });
    },
    { rootMargin: "-15% 0px -65% 0px" },
  );
  ["evidence", "studies", "roadmap"].forEach((id) => observer.observe(document.getElementById(id)));
  loadSnapshot();
})();
