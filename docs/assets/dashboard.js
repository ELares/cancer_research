/* Optional legacy workspace. Python, CDN assets and data load only on launch. */
(() => {
  "use strict";
  const RAW = "https://raw.githubusercontent.com/ELares/cancer_research/main";
  const STLITE = "https://cdn.jsdelivr.net/npm/@stlite/mountable@0.67.0/build";
  const FILES = [
    "scripts/dashboard.py",
    "scripts/dashboard_data.py",
    "corpus/INDEX.jsonl",
    "analysis/census-diagnostic-chains.json",
    "analysis/census-evidence-design.json",
    "analysis/census-mechanism-growth.json",
    "analysis/census-mechanism-profile.json",
    "analysis/census-mechanism-sites.json",
    "analysis/uncertainty-intervals-report.md",
  ];
  const READY_TEXT = "All dashboard panels are ready.";
  const button = document.getElementById("launch-dashboard");
  const status = document.getElementById("dashboard-status");
  const container = document.getElementById("dashboard-container");
  let attempt = 0;
  let observer;
  let timer;
  let downloads;

  function validateFile(file, content) {
    if (!content.trim()) throw new Error("The downloaded file is empty.");
    const isObject = (value) => value !== null && typeof value === "object" && !Array.isArray(value);
    if (file.endsWith(".json")) {
      const data = JSON.parse(content);
      if (!isObject(data)) throw new Error("Expected a JSON object.");
    } else if (file.endsWith(".jsonl")) {
      for (const line of content
        .trim()
        .split(/\r?\n/)
        .filter((value) => value.trim())) {
        const record = JSON.parse(line);
        if (!isObject(record)) throw new Error("Expected an archive record object.");
      }
    }
  }

  function setStatus(message, state) {
    status.textContent = message;
    status.dataset.state = state;
  }

  async function launchDashboard() {
    const currentAttempt = ++attempt;
    observer?.disconnect();
    clearTimeout(timer);
    downloads?.abort();
    const controller = new AbortController();
    downloads = controller;
    button.disabled = true;
    button.textContent = "Launching dashboard…";
    setStatus(
      "Loading Python and the committed data. The dashboard is not ready yet; first load can take 30–60 seconds.",
      "loading",
    );
    container.hidden = false;
    const frame = document.createElement("iframe");
    frame.title = "Research dashboard: census summaries, historical archive and read-only simulation report";
    container.replaceChildren(frame);
    const doc = frame.contentDocument;
    doc.open();
    doc.write(
      '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Research dashboard</title><style>html,body{margin:0;height:100%;font-family:system-ui,sans-serif}#stlite-root{height:100vh}.loading{padding:2rem;color:#52665a}</style></head><body><div id="stlite-root"><p class="loading">Preparing the research dashboard…</p></div></body></html>',
    );
    doc.close();
    let failed = false;

    function fail(message) {
      if (currentAttempt !== attempt) return;
      failed = true;
      clearTimeout(timer);
      controller.abort();
      button.disabled = false;
      button.textContent = "Retry dashboard";
      setStatus(message + " You can retry or return to the atlas above.", "error");
    }

    frame.contentWindow.addEventListener("error", (event) => {
      if (event.message) fail("The dashboard encountered an error: " + event.message);
    });
    frame.contentWindow.addEventListener("unhandledrejection", (event) => {
      fail("The dashboard could not finish loading: " + (event.reason?.message || String(event.reason)));
    });
    observer = new MutationObserver(() => {
      if (currentAttempt !== attempt) return;
      // Pinned stlite reports Pyodide/worker boot failures in an error toast,
      // without reliably forwarding them to the iframe's window events.
      const bootError = doc.querySelector(".Toastify__toast--error");
      if (bootError?.textContent.includes("Error during booting up")) {
        fail("The dashboard could not finish loading Python. Details are shown below.");
        return;
      }
      const exception = doc.querySelector('[data-testid="stException"]');
      if (exception) {
        fail("The dashboard reported an error. Details are shown below.");
        return;
      }
      // This marker is written by Python after all three panels finish. A
      // Streamlit shell or skeleton alone is not a successfully loaded app.
      if (!failed && doc.body.textContent.includes(READY_TEXT)) {
        clearTimeout(timer);
        button.disabled = true;
        button.textContent = "Dashboard launched";
        setStatus("Dashboard ready. Article filters apply to the 4,830-record historical archive.", "ready");
      }
    });
    observer.observe(doc.body, {
      childList: true,
      subtree: true,
      characterData: true,
    });
    timer = setTimeout(() => {
      if (currentAttempt !== attempt || status.dataset.state !== "loading") return;
      button.disabled = false;
      button.textContent = "Retry dashboard";
      setStatus(
        "The dashboard has not finished loading. Downloads may be slow or blocked; you can retry or keep waiting.",
        "error",
      );
    }, 90000);

    try {
      // stlite's worker does not propagate every file-download failure to this
      // window. Fetch here so HTTP/network errors immediately enable Retry, and
      // reject invalid data before Python can silently omit an affected panel.
      const files = Object.fromEntries(
        await Promise.all(
          FILES.map(async (file) => {
            try {
              const response = await fetch(`${RAW}/${file}`, {
                signal: controller.signal,
              });
              if (!response.ok) throw new Error(`HTTP ${response.status}`);
              const content = await response.text();
              validateFile(file, content);
              return [file, content];
            } catch (error) {
              throw new Error(`Could not load dashboard file ${file}: ${error.message || String(error)}`);
            }
          }),
        ),
      );
      if (currentAttempt !== attempt) return;
      await new Promise((resolve, reject) => {
        const css = doc.createElement("link");
        css.rel = "stylesheet";
        css.href = `${STLITE}/stlite.css`;
        css.onload = resolve;
        css.onerror = () => reject(new Error("The dashboard stylesheet could not be downloaded."));
        doc.head.appendChild(css);
      });
      await new Promise((resolve, reject) => {
        const script = doc.createElement("script");
        script.src = `${STLITE}/stlite.js`;
        script.onload = resolve;
        script.onerror = () => reject(new Error("The Python dashboard loader could not be downloaded."));
        doc.head.appendChild(script);
      });
      if (currentAttempt !== attempt) return;
      await frame.contentWindow.stlite.mount(
        {
          requirements: ["pandas"],
          entrypoint: "scripts/dashboard.py",
          files,
        },
        doc.getElementById("stlite-root"),
      );
    } catch (error) {
      fail(error.message || "The dashboard could not start.");
    }
  }

  button.addEventListener("click", launchDashboard);
})();
