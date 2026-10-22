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
# macOS `mktemp -d` with no template ignores TMPDIR and uses the per-user Darwin
# temp dir, which is unwritable under some sandboxes. Honour TMPDIR explicitly.
TMP=$(mktemp -d "${TMPDIR:-/tmp}/seqeval-smoke.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

fail() { echo "FAIL: $1" >&2; exit 1; }
pass() { echo "  ok: $1"; }

command -v uv >/dev/null 2>&1 || fail "uv not found on PATH"
[ -f "$FASTA" ] || fail "missing fixture $FASTA"

VL=DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK
VH=EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS
HC_LALA=$(awk '/HC_LALA/{f=1;next} /^>/{f=0} f' "$FASTA" | tr -d '\n')
HC_WT=$(awk '/HC_full/{f=1;next} /^>/{f=0} f' "$FASTA" | tr -d '\n')
HC_IGG4=$(awk '/IgG4_WT/{f=1;next} /^>/{f=0} f' "$FASTA" | tr -d '\n')
HC_IGG4_S228P=$(awk '/IgG4_S228P/{f=1;next} /^>/{f=0} f' "$FASTA" | tr -d '\n')

echo "[1/6] analyze_properties.py"
uv run scripts/analyze_properties.py "$FASTA" --json > "$TMP/props.json"
python3 -c "import json,sys; d=json.load(open('$TMP/props.json')); assert len(d['chains'])==7, d; assert d['chains'][0]['a280_oxidized']>0" \
  || fail "properties json malformed"
pass "7 chains, A280 computed"

echo "[2/6] find_cdrs.py (VL)"
uv run scripts/find_cdrs.py --vl "$VL" --name VL --json > "$TMP/vl.json"
python3 -c "import json; d=json.load(open('$TMP/vl.json')); assert d['cdr3']=='QQHYTTPPT', d; assert d['domain_type']=='VL'" \
  || fail "VL CDR3 wrong"
pass "VL CDRs correct (QQHYTTPPT)"

echo "[3/6] find_mutations.py (WT + LALA annotation)"
uv run scripts/find_mutations.py "$HC_WT" --label WT --json > "$TMP/wt.json"
python3 -c "import json; d=json.load(open('$TMP/wt.json')); assert d['substitutions']==[], d" \
  || fail "WT IgG1 should have 0 substitutions"
python3 -c "import json; d=json.load(open('$TMP/wt.json')); assert d['isotype']=='IgG1', d['isotype']" \
  || fail "WT not detected as IgG1"
uv run scripts/find_mutations.py "$HC_LALA" --label LALA --json > "$TMP/lala.json"
python3 -c "import json; d=json.load(open('$TMP/lala.json')); m={s['mutation'] for s in d['substitutions']}; assert m=={'L234A','L235A'}, m; assert 'LALA' in d['variants_present']" \
  || fail "LALA not detected/annotated"
pass "WT=0 subs (IgG1); LALA detected and annotated"

echo "[3b/6] find_mutations.py isotype awareness (IgG4 + S228P)"
uv run scripts/find_mutations.py "$HC_IGG4" --label IgG4 --json > "$TMP/igg4.json"
python3 -c "import json; d=json.load(open('$TMP/igg4.json')); assert d['isotype']=='IgG4', d['isotype']; assert d['substitutions']==[], d['substitutions']" \
  || fail "IgG4 WT not auto-detected as clean IgG4"
uv run scripts/find_mutations.py "$HC_IGG4_S228P" --label IgG4P --json > "$TMP/igg4p.json"
python3 -c "import json; d=json.load(open('$TMP/igg4p.json')); assert d['isotype']=='IgG4'; subs=d['substitutions']; assert len(subs)==1 and subs[0]['mutation']=='S228P' and subs[0]['eu']==228, subs" \
  || fail "IgG4 S228P not detected at EU 228"
pass "IgG4 auto-detected; S228P called at EU 228"

echo "[4/6] scan_liabilities.py"
uv run scripts/scan_liabilities.py "$FASTA" --json > "$TMP/liab.json"
python3 -c "import json; d=json.load(open('$TMP/liab.json')); assert 'Trastuzumab_VH' in d; assert 'findings' in d['Trastuzumab_VH']; assert not d['Trastuzumab_VH']['cdr_annotated']; assert all(not f['escalated'] for f in d['Trastuzumab_VH']['findings'])" \
  || fail "liabilities json malformed or escalated without a CDR map"
pass "liabilities scanned for all chains; nothing escalated without --cdr-json"

echo "[4b/6] build_cdr_map.py + CDR-escalated liabilities"
uv run scripts/find_cdrs.py --vh "$VH" --name Trastuzumab_VH --json > "$TMP/vh.json"
uv run scripts/build_cdr_map.py Trastuzumab_VH="$TMP/vh.json" > "$TMP/cdrmap.json"
python3 -c "import json; m=json.load(open('$TMP/cdrmap.json')); e=m['Trastuzumab_VH']; assert [x['name'] for x in e]==['CDR-H1','CDR-H2','CDR-H3'], e" \
  || fail "build_cdr_map did not emit CDR-H1/2/3"
uv run scripts/scan_liabilities.py --seq "$VH" --name Trastuzumab_VH --cdr-json "$TMP/cdrmap.json" --json > "$TMP/liab_cdr.json"
python3 -c "
import json
d=json.load(open('$TMP/liab_cdr.json'))['Trastuzumab_VH']
assert d['cdr_annotated'], d
esc=[f for f in d['findings'] if f['escalated']]
assert esc, 'no finding escalated despite a CDR map'
assert all(f['in_cdr'] and f['cdr'] for f in esc), esc
dg=[f for f in d['findings'] if f['motif']=='DG' and f.get('cdr')=='CDR-H3']
assert dg and dg[0]['severity']=='high' and dg[0]['base_severity']=='medium', dg
" || fail "CDR escalation did not fire"
pass "CDR map built; CDR-H3 isomerization escalated medium->high"
# conserved Fc N297 sequon must be labelled 'expected', not a liability
python3 -c "
import json; d=json.load(open('$TMP/liab.json'))
hc=d['Trastuzumab_HC_full_IgG1']['findings']
n297=[f for f in hc if f.get('fc_n297')]
assert len(n297)==1, n297
assert n297[0]['severity']=='expected', n297[0]
assert not n297[0].get('escalated'), 'N297 must never be escalated'
" || fail "Fc N297 not labelled expected"
pass "conserved Fc N297 reported as expected, exempt from escalation"

echo "[2b/6] find_cdrs.py CDR3 anchors (VYFCAR FR3 + Kabat H3 boundary)"
# Trastuzumab VH with FR3 mutated YYCSR->YFCAR: reproduces the real-world
# 'VYFCAR' case that the old literal-YYC anchor could not parse.
VH_YFC=EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYFCARWGGDGFYAMDYWGQGTLVTVSS
uv run scripts/find_cdrs.py --vh "$VH_YFC" --name YFC --json > "$TMP/yfc.json"
python3 -c "import json; d=json.load(open('$TMP/yfc.json')); assert d['cdr3']=='WGGDGFYAMDY', d" \
  || fail "VYFCAR-style FR3 not handled"
TVH=EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS
uv run scripts/find_cdrs.py --vh "$TVH" --name TVH --json > "$TMP/tvh.json"
python3 -c "import json; d=json.load(open('$TMP/tvh.json')); assert d['cdr3']=='WGGDGFYAMDY', d" \
  || fail "trastuzumab Kabat CDR-H3 boundary wrong"
pass "CDR-H3 anchored on Cys92; YFC FR3 works; Kabat boundary exact"
echo "[5/6] analyze_interfaces.py"
uv run scripts/analyze_interfaces.py examples/mini_complex.cif --antibody A --antigen C \
  --cdr-json examples/mini_cdrs.json --json > "$TMP/iface.json"
python3 -c "import json; d=json.load(open('$TMP/iface.json')); assert d['total_atom_contacts']>0; assert any(p['region']=='CDR-H3' for p in d['paratope'])" \
  || fail "interface analysis / CDR mapping failed"
pass "contacts found; paratope mapped to CDR-H3"

echo "[6/7] scan_surface_patches.py"
uv run scripts/scan_surface_patches.py examples/mini_patches.cif --chains A \
  --cdr-json examples/mini_patches_cdrs.json --json > "$TMP/patches.json"
python3 -c "
import json
d=json.load(open('$TMP/patches.json'))
pos=[p for p in d['patches'] if p['kind']=='positive_charge']
hyd=[p for p in d['patches'] if p['kind']=='hydrophobic']
assert d['method']=='surface-patch-v1' and d['not_tap'] and d['not_sap']
assert len(pos)==1 and pos[0]['severity']=='high' and pos[0]['n_residues']==3 and pos[0]['cdr_residues']==3, pos
assert len(hyd)==1 and hyd[0]['severity']=='medium' and hyd[0]['n_residues']==3 and hyd[0]['cdr_residues']==0, hyd
assert not d['coarse_model']
assert d['chains_in_file']==['A'] and d['chains_excluded_from_sasa']==[], d
" || fail "surface patch scan failed"
pass "CDR positive patch high; framework hydrophobic patch medium"

echo "[6b/7] surface patches exclude non-selected chains from SASA"
uv run scripts/scan_surface_patches.py examples/mini_complex.cif --chains A --json > "$TMP/patches_strip.json"
python3 -c "
import json
d=json.load(open('$TMP/patches_strip.json'))
assert d['chains']==['A'], d['chains']
assert 'C' in d['chains_in_file'], d['chains_in_file']
assert d['chains_excluded_from_sasa']==['C'], d['chains_excluded_from_sasa']
assert any('excluded before the SASA' in w for w in d['warnings']), d['warnings']
" || fail "antigen chain was not excluded from the SASA calculation"
pass "antigen chain C detached before SASA; exclusion recorded"

echo "[7/7] compile_report.py"
uv run scripts/compile_report.py --title "Smoke Test" \
  --properties "$TMP/props.json" --cdrs "$TMP/vl.json" \
  --mutations "$TMP/lala.json" --liabilities "$TMP/liab.json" \
  --patches "$TMP/patches.json" > "$TMP/report.md"
for phase in "Phase 1 — Molecule Triage" "Phase 2 — Chain Properties" \
             "Phase 3 — Domain Annotation" "Phase 4 — Engineering Audit" \
             "Phase 5 — Sequence Liabilities" "Phase 6 — Structural Developability" \
             "Phase 7 — Lineage" "Phase 8 — Verdict"; do
  grep -q "## $phase" "$TMP/report.md" || fail "report missing $phase"
done
grep -q "Binder Assessment" "$TMP/report.md" && fail "binder phase should have been removed"
grep -q "## Phase 9" "$TMP/report.md" && fail "report should end at Phase 8"
grep -q "LALA" "$TMP/report.md" || fail "report missing LALA annotation"
grep -q "PPC/PNC/PSH" "$TMP/report.md" || fail "report missing TAP caveat"
grep -q "## Provenance" "$TMP/report.md" || fail "report missing provenance block"
# Every phase must carry an Assessment or an explicit not-assessed line.
python3 -c "
import re,sys
txt=open('$TMP/report.md').read()
blocks=re.split(r'^## Phase ', txt, flags=re.M)[1:]
for b in blocks:
    name=b.splitlines()[0].strip()
    if '**Assessment**' not in b and 'Not assessed' not in b:
        sys.exit('phase without an assessment: '+name)
" || fail "a phase emitted findings with no assessment"
grep -qE "STOP|ENGINEER|PROCEED" "$TMP/report.md" || fail "verdict has no overall call"
pass "8 phases present, each with an assessment, plus verdict and provenance"

echo "[7b/7] verdict escalates a construct error to STOP"
python3 -c "
import json
json.dump({'label':'ArmA','isotype':'IgG1','substitutions':[{'eu':366,'domain':'CH3','mutation':'T366W','note':'KiH knob'}]},open('$TMP/knob.json','w'))
json.dump({'label':'ArmB','isotype':'IgG1','substitutions':[{'eu':366,'domain':'CH3','mutation':'T366W','note':'KiH knob'}]},open('$TMP/knob2.json','w'))
"
uv run scripts/compile_report.py --title "KiH check" --mutations "$TMP/knob.json" "$TMP/knob2.json" > "$TMP/kih.md"
grep -q "Construct error" "$TMP/kih.md" || fail "two knobs not reported as a construct error"
grep -q "STOP" "$TMP/kih.md" || fail "two knobs did not produce a STOP verdict"
pass "knob-on-both-arms reported as construct error and STOP"

echo "[8/8] select_structure_model.py routing"
printf "%s" "$VH" > "$TMP/vh.seq"
printf "%s" "$VL" > "$TMP/vl.seq"
uv run scripts/select_structure_model.py --cdrs "$TMP/vh.json" "$TMP/vl.json" \
  --seqs VH="$TMP/vh.seq" VL="$TMP/vl.seq" --out-fasta "$TMP/fv.fasta" --json > "$TMP/route.json"
python3 -c "
import json
d=json.load(open('$TMP/route.json'))
assert d['modality']=='FAB', d
assert d['model']=='Chai-1' and d['skill']=='chai', d
assert d['gpu_required'] and d['job_class']=='cofold_heavy', d
assert [c['chain'] for c in d['chains']]==['A','B'], d['chains']
assert d['antibody_chains']==['A','B'] and d['antigen_chains']==[], d
" || fail "Fab did not route to Chai-1"
VHH_SEQ=QVQLVESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYFCAAKYGSSWYLDAWGQGTQVTVSS
uv run scripts/find_cdrs.py --vhh "$VHH_SEQ" --name Nb --json > "$TMP/vhh.json"
printf "%s" "$VHH_SEQ" > "$TMP/vhh.seq"
uv run scripts/select_structure_model.py --cdrs "$TMP/vhh.json" \
  --seqs VHH="$TMP/vhh.seq" --out-fasta "$TMP/nb.fasta" --json > "$TMP/route_nb.json"
python3 -c "
import json
d=json.load(open('$TMP/route_nb.json'))
assert d['modality']=='VHH', d
assert d['model']=='Boltz-1' and d['skill']=='boltz', d
assert d['antibody_chains']==['A'], d
" || fail "VHH did not route to Boltz-1"
grep -q "^>VHH|A" "$TMP/nb.fasta" || fail "staged FASTA missing chain-id header"
pass "Fab -> Chai-1, VHH -> Boltz-1; chain IDs staged"

echo "[preflight] preflight.py"
uv run scripts/preflight.py --json > "$TMP/pre.json" || true
python3 -c "
import json; d=json.load(open('$TMP/pre.json'))
r=d['recommendations']
assert set(r)=={'cpu_light','ab_structure','plm_embed','cofold_heavy'}, r
assert r['cpu_light']['backend']=='local'
assert d['tooling']['uv'], 'uv not detected'
mods={p['module'] for p in d['packages']}
assert {'Bio','numpy','modal'} <= mods, mods
assert d['required_ok'] is True, d['action_items']
" || fail "preflight output malformed"
pass "preflight reports tooling, packages, GPU and backends"

echo "[2d/6] CDR3 must not truncate when CDR3 itself contains a cysteine"
TVH_C=EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDCGFYAMDYWGQGTLVTVSS
TVL_C=DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYCTPPTFGQGTKVEIK
uv run scripts/find_cdrs.py --vh "$TVH_C" --name hcys --json > "$TMP/hcys.json"
python3 -c "import json; d=json.load(open('$TMP/hcys.json')); assert d['cdr3']=='WGGDCGFYAMDY', d" \
  || fail "CDR-H3 truncated by an internal cysteine"
uv run scripts/find_cdrs.py --vl "$TVL_C" --name lcys --json > "$TMP/lcys.json"
python3 -c "import json; d=json.load(open('$TMP/lcys.json')); assert d['cdr3']=='QQHYCTPPT', d" \
  || fail "CDR-L3 truncated by an internal cysteine"
# a CDR3 containing the J motif internally must not end the loop early
TVH_J=EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGQGTYAMDYWGQGTLVTVSS
uv run scripts/find_cdrs.py --vh "$TVH_J" --name hj --json > "$TMP/hj.json"
python3 -c "import json; d=json.load(open('$TMP/hj.json')); assert d['cdr3']=='WGQGTYAMDY', d" \
  || fail "internal J motif truncated CDR-H3"
pass "forward Cys anchor + last J match survive Cys-in-CDR3 and internal WGQGT"

echo "[2c/6] ANARCI backend-missing fallback"
python3 -I tests/test_anarci_fallback.py >/dev/null \
  || fail "find_cdrs does not fall back when ANARCI backend is missing"
pass "abnumber-without-anarci degrades to regex anchors with guidance"

echo "[modal] cloud GPU app definitions"
python3 -m py_compile modal/modal_chai1.py modal/modal_boltz.py \
  || fail "modal app syntax error"
grep -q "local_entrypoint" modal/modal_chai1.py && grep -q "local_entrypoint" modal/modal_boltz.py \
  || fail "modal apps missing local_entrypoint"
pass "modal Chai-1/Boltz apps compile and expose entrypoints"

echo "[pre-flight b] blocked cofold_heavy emits a remediation ladder"
python3 -c "
import json
d=json.load(open('$TMP/pre.json'))
if 'cofold_heavy' in d['blocked_job_classes']:
    fixes=d['remediation']
    assert fixes, 'blocked but no remediation offered'
    assert fixes[0]['option']=='Modal' and fixes[0]['recommended'], fixes[0]
    assert not fixes[0]['changes_model'], fixes[0]
    assert any(f['changes_model'] for f in fixes), 'model-changing options not flagged'
else:
    assert d['remediation']==[], 'remediation offered for a machine that is not blocked'
" || fail "remediation ladder wrong"
uv run scripts/preflight.py --json --modal yes > "$TMP/res_modal.json"
python3 -c "
import json
d=json.load(open('$TMP/res_modal.json'))
assert d['modal_available'] and d['modal_source']=='asserted', d
assert d['recommendations']['cofold_heavy']['backend']=='modal', d['recommendations']['cofold_heavy']
assert d['blocked_job_classes']==[] and d['remediation']==[], d
" || fail "--modal yes did not route cofold_heavy to modal"
pass "ladder offered when blocked; --modal yes routes cofold_heavy to Modal"

echo "ALL SMOKE TESTS PASSED"
