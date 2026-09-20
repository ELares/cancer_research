/* Exercise the shipped controller without CDN, Python or browser dependencies. */
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const test = require("node:test");
const vm = require("node:vm");

const root = path.resolve(__dirname, "..");
const source = fs.readFileSync(path.join(root, "docs/assets/dashboard.js"), "utf8");
const raw = "https://raw.githubusercontent.com/ELares/cancer_research/main/";

function harness(fetchOverride) {
  const requests = [];
  const mounts = [];
  const frames = [];
  const timers = [];
  const observers = [];
  const button = {
    addEventListener: (_event, callback) => {
      button.launch = callback;
    },
  };
  const status = { dataset: {} };
  const container = { replaceChildren: (frame) => frames.push(frame) };
  const elements = {
    "launch-dashboard": button,
    "dashboard-status": status,
    "dashboard-container": container,
  };
  const document = {
    getElementById: (id) => elements[id],
    createElement: (tag) => {
      assert.equal(tag, "iframe");
      const body = { textContent: "" };
      const doc = {
        body,
        nodes: new Map(),
        open() {},
        write() {},
        close() {},
        createElement: (type) => ({ tagName: type }),
        head: {
          appendChild: (element) => queueMicrotask(() => element.onload()),
        },
        getElementById: () => ({}),
        querySelector: (selector) => doc.nodes.get(selector) ?? null,
      };
      return {
        contentDocument: doc,
        contentWindow: {
          addEventListener() {},
          stlite: {
            mount: async (options) => {
              mounts.push(options);
            },
          },
        },
      };
    },
  };
  vm.runInNewContext(source, {
    document,
    AbortController,
    MutationObserver: class {
      constructor(callback) {
        this.callback = callback;
        observers.push(this);
      }
      observe() {}
      disconnect() {}
    },
    setTimeout: (callback) => {
      timers.push(callback);
      return callback;
    },
    clearTimeout() {},
    fetch: async (url, options) => {
      const file = url.slice(raw.length);
      requests.push({ file, signal: options.signal });
      const result = await fetchOverride?.(file, options, requests);
      return (
        result ?? {
          ok: true,
          text: async () => fs.readFileSync(path.join(root, file), "utf8"),
        }
      );
    },
  });
  return { button, status, requests, mounts, frames, timers, observers };
}

test("files download only after launch and mount as checked text contents", async () => {
  const app = harness();
  assert.equal(app.requests.length, 0);
  assert.equal(app.frames.length, 0);
  await app.button.launch();
  assert.equal(app.mounts.length, 1);
  const files = app.mounts[0].files;
  assert.equal(Object.keys(files).length, 9);
  for (const [file, content] of Object.entries(files)) {
    assert.equal(content, fs.readFileSync(path.join(root, file), "utf8"));
  }
  assert.equal(app.status.dataset.state, "loading");
  app.frames[0].contentDocument.body.textContent = "All dashboard panels are ready.";
  app.observers[0].callback();
  assert.equal(app.status.dataset.state, "ready");
  assert.equal(app.button.disabled, true);
});

test("an HTTP failure enables immediate retry and retry fetches every file again", async () => {
  let fail = true;
  const app = harness((file) => {
    if (fail && file === "analysis/census-evidence-design.json") {
      return { ok: false, status: 404, text: async () => "Not Found" };
    }
  });
  await app.button.launch();
  assert.equal(app.status.dataset.state, "error");
  assert.match(app.status.textContent, /census-evidence-design.json: HTTP 404/);
  assert.equal(app.button.disabled, false);
  assert.equal(app.button.textContent, "Retry dashboard");
  assert.equal(app.mounts.length, 0);
  assert.ok(app.requests.every((request) => request.signal.aborted));
  fail = false;
  await app.button.launch();
  assert.equal(app.mounts.length, 1);
  assert.equal(app.frames.length, 2);
  assert.equal(app.requests.length, 18);
  assert.equal(app.status.dataset.state, "loading");
});

test("a rejected request reports its file and does not mount an incomplete app", async () => {
  const app = harness((file) => {
    if (file === "corpus/INDEX.jsonl") throw new TypeError("Failed to fetch");
  });
  await app.button.launch();
  assert.equal(app.status.dataset.state, "error");
  assert.match(app.status.textContent, /corpus\/INDEX.jsonl: Failed to fetch/);
  assert.equal(app.button.disabled, false);
  assert.equal(app.mounts.length, 0);
});

test("a worker bootstrap error immediately enables retry and cannot become ready", async () => {
  const app = harness();
  await app.button.launch();
  const doc = app.frames[0].contentDocument;
  doc.nodes.set(".Toastify__toast--error", {
    textContent: "Error during booting up\nFailed to fetch",
  });
  app.observers[0].callback();
  assert.equal(app.status.dataset.state, "error");
  assert.match(app.status.textContent, /could not finish loading Python/);
  assert.equal(app.button.disabled, false);
  doc.nodes.clear();
  doc.body.textContent = "All dashboard panels are ready.";
  app.observers[0].callback();
  assert.equal(app.status.dataset.state, "error");
  await app.button.launch();
  app.frames[1].contentDocument.body.textContent = "All dashboard panels are ready.";
  app.observers[1].callback();
  assert.equal(app.status.dataset.state, "ready");
});

test("invalid or empty downloaded data cannot produce a ready app", async (t) => {
  for (const [file, content] of [
    ["analysis/census-evidence-design.json", "not JSON"],
    ["analysis/census-evidence-design.json", "null"],
    ["corpus/INDEX.jsonl", '{"pmid":"1"}\nnot JSON'],
    ["corpus/INDEX.jsonl", "[]"],
    ["analysis/uncertainty-intervals-report.md", "  \n"],
  ]) {
    await t.test(`${file}: ${JSON.stringify(content)}`, async () => {
      const app = harness((requested) => (requested === file ? { ok: true, text: async () => content } : undefined));
      await app.button.launch();
      assert.equal(app.status.dataset.state, "error");
      assert.ok(app.status.textContent.includes(file));
      assert.equal(app.button.disabled, false);
      assert.equal(app.mounts.length, 0);
    });
  }
});

test("retry aborts a slow attempt and ignores its late completion", async () => {
  let release;
  let first = true;
  const app = harness((file) => {
    if (first && file === "scripts/dashboard.py") {
      first = false;
      return new Promise((resolve) => {
        release = resolve;
      });
    }
  });
  const slowAttempt = app.button.launch();
  app.timers[0]();
  assert.equal(app.button.disabled, false);
  assert.equal(app.status.dataset.state, "error");
  await app.button.launch();
  assert.equal(app.mounts.length, 1);
  assert.ok(app.requests.slice(0, 9).every((request) => request.signal.aborted));
  release({ ok: true, text: async () => "print('old attempt')" });
  await slowAttempt;
  assert.equal(app.mounts.length, 1);
  assert.notEqual(app.mounts[0].files["scripts/dashboard.py"], "print('old attempt')");
  assert.equal(app.status.dataset.state, "loading");
});

test("a silent bootstrap stall offers retry but a slow app may still become ready", async () => {
  const app = harness();
  await app.button.launch();
  app.timers[0]();
  assert.equal(app.status.dataset.state, "error");
  assert.match(app.status.textContent, /slow or blocked/);
  assert.equal(app.button.disabled, false);
  app.frames[0].contentDocument.body.textContent = "All dashboard panels are ready.";
  app.observers[0].callback();
  assert.equal(app.status.dataset.state, "ready");
  assert.equal(app.button.disabled, true);
});
