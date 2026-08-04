#!/usr/bin/env bash
# The atlas pipeline, in dependency order (#ATLAS).
#
# WHY THIS EXISTS
# ---------------
# The chain that built the census lived in a gitignored log directory and
# covered four steps: baseline, coverage, fulltext, relations, graph. Every
# quality layer added since -- ambiguity, disambiguation, domain-sense
# validation, impact, co-mention, and the audits -- had to be run by hand, in an
# order that is not obvious and that matters.
#
# It matters because the layers feed each other:
#
#   relations  rebuilds the ENTITY files, which invalidates the ambiguity scan
#   ambiguity  produces the blocklist that disambiguation and co-mention consume
#   disambiguate produces the per-paper corrections the GRAPH BUILD applies
#   graph      must therefore be rebuilt AFTER disambiguation, not before
#   co-mention consumes the blocklist and must be --rebuild after a corpus change
#
# Running these out of order does not fail loudly. It produces a graph built on
# stale corrections, or a co-mention table counting collided aliases, and both
# look fine.
#
# NOT a cron job. Several steps take hours and one needs a decision (see the
# NEWS section). Run it deliberately, or run the phases individually.
#
# Usage:
#   bash scripts/atlas_pipeline.sh census     # baseline -> relations (hours)
#   bash scripts/atlas_pipeline.sh quality    # ambiguity -> corrected graph
#   bash scripts/atlas_pipeline.sh mine       # the analyses that read the graph
#   bash scripts/atlas_pipeline.sh audit      # citation + news integrity checks
set -euo pipefail

PY="${PY:-.venv/bin/python}"
cd "$(dirname "$0")/.."
PHASE="${1:-}"

log() { printf '\n=== %s === %s\n' "$1" "$(date)"; }

census() {
  log "census: PubMed baseline -> cancer records"
  $PY scripts/atlas_baseline.py
  log "coverage: frozen corpus measured against the census"
  $PY scripts/atlas_coverage.py
  log "fulltext: PMC bulk (hours; recency ceiling documented in the script)"
  $PY scripts/atlas_fulltext.py
  log "relations: PubTator3 bulk -- REBUILDS THE ENTITY FILES"
  $PY scripts/atlas_relations.py
}

quality() {
  # Order is load-bearing from here down.
  log "ambiguity: measure sense collisions -> blocklist (~40 min, network)"
  $PY scripts/atlas_ambiguity.py --top 400
  log "domain sense: verify the curated senses against the corpus"
  $PY scripts/atlas_domain_sense.py --sample 800
  log "disambiguate: per-paper FSP1 senses -> corrections"
  $PY scripts/atlas_disambiguate.py
  log "graph: rebuild, APPLYING the corrections produced above"
  $PY scripts/atlas_graph.py --build
  log "impact: how far the collisions actually reach"
  $PY scripts/atlas_ambiguity_impact.py
  log "co-mention: full-text sentences (~2 h). --rebuild clears table AND manifest"
  $PY scripts/atlas_comention.py --rebuild
}

mine() {
  log "module support"      ; $PY scripts/atlas_module_support.py
  log "contradictions"      ; $PY scripts/atlas_contradictions.py
  log "contradiction quality"; $PY scripts/atlas_contradiction_quality.py
  log "emergence"           ; $PY scripts/atlas_emergence.py
  log "emergence error"     ; $PY scripts/atlas_emergence_error.py
  log "discovery evaluation"; $PY scripts/atlas_discovery_eval.py
  log "entity audit"        ; $PY scripts/atlas_entity_audit.py
}

audit() {
  log "citation audit: do the cited papers exist and match their claims"
  $PY scripts/atlas_citation_audit.py
  log "news verification audit: does the 'verified' label hold"
  $PY scripts/news_verification_audit.py
  # Deliberately NOT re-running scripts/verify_news_claims.py. That rewrites
  # claim statuses and their credibility scores (the post-fix re-run withdrew 30
  # of 44 verifications and moved every article's score), which is a data change
  # to review on its own rather than a pipeline side effect. It also depends on
  # PubMed's live index, so it is not reproducible the way the rest of this
  # pipeline is.
}

case "$PHASE" in
  census)  census ;;
  quality) quality ;;
  mine)    mine ;;
  audit)   audit ;;
  *) sed -n '/^# Usage:/,/^set -euo/p' "$0" | head -8; exit 2 ;;
esac
log "$PHASE complete"
