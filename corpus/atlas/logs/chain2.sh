#!/bin/bash
cd /Users/ezequiellares/research
L=corpus/atlas/logs
.venv/bin/python scripts/atlas_baseline.py >> $L/chain2.log 2>&1
echo "[chain2] census-v2 done $(date)" >> $L/chain2.log
.venv/bin/python scripts/atlas_coverage.py >> $L/chain2.log 2>&1
echo "[chain2] coverage done $(date)" >> $L/chain2.log
.venv/bin/python scripts/atlas_fulltext.py >> $L/chain2.log 2>&1
echo "[chain2] fulltext done $(date)" >> $L/chain2.log
.venv/bin/python scripts/atlas_relations.py >> $L/chain2.log 2>&1
echo "[chain2] relations done $(date)" >> $L/chain2.log
.venv/bin/python scripts/atlas_graph.py --build >> $L/chain2.log 2>&1
echo "[chain2] graph rebuilt $(date)" >> $L/chain2.log
