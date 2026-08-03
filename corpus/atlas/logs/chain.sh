#!/bin/bash
cd /Users/ezequiellares/research
L=corpus/atlas/logs
while pgrep -f "atlas_baseline.py" >/dev/null 2>&1; do sleep 60; done
echo "[chain] census done $(date)" >> $L/chain.log
.venv/bin/python scripts/atlas_coverage.py  >> $L/chain.log 2>&1
echo "[chain] coverage done $(date)" >> $L/chain.log
.venv/bin/python scripts/atlas_fulltext.py  >> $L/chain.log 2>&1
echo "[chain] fulltext done $(date)" >> $L/chain.log
.venv/bin/python scripts/atlas_relations.py >> $L/chain.log 2>&1
echo "[chain] relations done $(date)" >> $L/chain.log
