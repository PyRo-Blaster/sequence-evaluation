#!/usr/bin/env bash
# Smoke test for the sequence-evaluation scripts.
# Runs every tool against the committed trastuzumab fixtures and checks key
# outputs. Intended for CI and for a SessionStart hook on Claude Code on the web.
#
#   bash tests/smoke_test.sh
#
# Requires `uv` on PATH (https://docs.astral.sh/uv/). Dependencies are pulled in
# automatically via each script's inline metadata.
set -euo pipefail

cd "$(dirname "$0")/.."
FASTA=examples/trastuzumab.fasta
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

fail() { echo "FAIL: $1" >&2; exit 1; }
pass() { echo "  ok: $1"; }

command -v uv >/dev/null 2>&1 || fail "uv not found on PATH"
[ -f "$FASTA" ] || fail "missing fixture $FASTA"

VL=DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK
HC_LALA=$(awk '/HC_LALA/{f=1;next} /^>/{f=0} f' "$FASTA" | tr -d '\n')
HC_WT=$(awk '/HC_full/{f=1;next} /^>/{f=0} f' "$FASTA" | tr -d '\n')

echo "[1/6] analyze_properties.py"
uv run scripts/analyze_properties.py "$FASTA" --json > "$TMP/props.json"
python3 -c "import json,sys; d=json.load(open('$TMP/props.json')); assert len(d['chains'])==5, d; assert d['chains'][0]['a280_oxidized']>0" \
  || fail "properties json malformed"
pass "5 chains, A280 computed"

echo "[2/6] find_cdrs.py (VL)"
uv run scripts/find_cdrs.py --vl "$VL" --name VL --json > "$TMP/vl.json"
python3 -c "import json; d=json.load(open('$TMP/vl.json')); assert d['cdr3']=='QQHYTTPPT', d; assert d['domain_type']=='VL'" \
  || fail "VL CDR3 wrong"
pass "VL CDRs correct (QQHYTTPPT)"

echo "[3/6] find_mutations.py (WT + LALA annotation)"
uv run scripts/find_mutations.py "$HC_WT" --label WT --json > "$TMP/wt.json"
python3 -c "import json; d=json.load(open('$TMP/wt.json')); assert d['substitutions']==[], d" \
  || fail "WT IgG1 should have 0 substitutions"
uv run scripts/find_mutations.py "$HC_LALA" --label LALA --json > "$TMP/lala.json"
python3 -c "import json; d=json.load(open('$TMP/lala.json')); m={s['mutation'] for s in d['substitutions']}; assert m=={'L234A','L235A'}, m; assert 'LALA' in d['variants_present']" \
  || fail "LALA not detected/annotated"
pass "WT=0 subs; LALA detected and annotated"

echo "[4/6] scan_liabilities.py"
uv run scripts/scan_liabilities.py "$FASTA" --json > "$TMP/liab.json"
python3 -c "import json; d=json.load(open('$TMP/liab.json')); assert 'Trastuzumab_VH' in d; assert 'findings' in d['Trastuzumab_VH']" \
  || fail "liabilities json malformed"
pass "liabilities scanned for all chains"

echo "[5/6] analyze_interfaces.py"
uv run scripts/analyze_interfaces.py examples/mini_complex.cif --antibody A --antigen C \
  --cdr-json examples/mini_cdrs.json --json > "$TMP/iface.json"
python3 -c "import json; d=json.load(open('$TMP/iface.json')); assert d['total_atom_contacts']>0; assert any(p['region']=='CDR-H3' for p in d['paratope'])" \
  || fail "interface analysis / CDR mapping failed"
pass "contacts found; paratope mapped to CDR-H3"

echo "[6/6] compile_report.py"
uv run scripts/compile_report.py --title "Smoke Test" \
  --properties "$TMP/props.json" --cdrs "$TMP/vl.json" \
  --mutations "$TMP/lala.json" --liabilities "$TMP/liab.json" \
  --interfaces "$TMP/iface.json" > "$TMP/report.md"
grep -q "## Chain Properties" "$TMP/report.md" || fail "report missing properties section"
grep -q "LALA" "$TMP/report.md" || fail "report missing LALA annotation"
pass "report compiled with all sections"

echo "ALL SMOKE TESTS PASSED"
