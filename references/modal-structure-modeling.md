# Running structure modelling on Modal

How a `cofold_heavy` job reaches a GPU when this machine has none. Read
`SKILL.md` → Stage 4 first for *which* model runs; this file is only about
*where* it runs.

> **Status: not verified end to end.** No Modal target was enabled when this was
> written, so the image/volume/egress facts below are read from the shipped
> `remote-compute-modal` env definitions, and the submit shapes follow that
> skill's documented API — but no job has been dispatched through this path yet.
> Treat the first run as a smoke test.

## The escalation chain

```
preflight.py                   (local first)
   cofold_heavy -> local       GPU with >=16 GB VRAM present       -> run here
   cofold_heavy -> modal       no local GPU, Modal enabled         -> this file
   cofold_heavy -> web-api     no GPU, no Modal, network available -> hosted server
   cofold_heavy -> blocked     none of the above                   -> ladder printed
```

The pre-flight cannot see a host-managed Modal target — that one has no local
CLI — so after enabling it, re-run with `--modal yes`:

```bash
uv run scripts/preflight.py --modal yes
```

`cofold_heavy` then routes to `modal` instead of `blocked`. Enabling Modal is
the user's decision (it bills to their account); the pre-flight prints the ladder
and stops rather than configuring anything itself.

## Enabling the target

*Settings → Compute → Modal.* An agent cannot enable it:
`host.compute.create('modal')` raises `host.compute.NotFound` listing the
targets that do exist. Confirm with `list_compute` once enabled.

## Which bundled env carries which model

Both predictors already ship as `remote-compute-modal` envs — build the image
rather than hand-rolling one.

| Model | Env | Pin | GPU default | Weights volume |
|-------|-----|-----|-------------|----------------|
| Chai-1 (Fab / Fv) | `proteomics_gpu` | `chai_lab==0.6.1`, torch 2.5.1 | `A100-80GB` | `/weights/chai` |
| Boltz (VHH) | `proteomics_boltz_gpu` | `boltz[cuda]==2.2.1` | `A100-80GB` | `/root/.boltz` |

**The shipped Boltz env is Boltz-2, not Boltz-1.** If the project specifically
requires Boltz-1, that needs an ad-hoc env with a `boltz==1.x` pin
(`build_env(path=…)`, and record the `egress:` ledger line yourself) — the
bundled path gives you Boltz-2.

Job-time egress each env declares — the sandbox is fenced to these plus the
provider policy:

- `proteomics_gpu`: `api.colabfold.com`, `huggingface.co`, `*.hf.co`, `chaiassets.com`
- `proteomics_boltz_gpu`: `api.colabfold.com`, `model-gateway.boltz.bio`, `huggingface.co`, `*.hf.co`

Weight downloads belong in the volume hydrate, not the job: `build_env(name,
hydrate=True)` populates the volume from the `compute_provider` kernel, which is
outside the fence.

## Two surfaces

Image building and job submission are different kernels on purpose:

- `compute_provider` tool → build the image, hydrate the weights volume. Returns
  an `im-…` reference. Capped at 30 min, no GPU.
- `repl` tool → `host.compute.create('modal', provider_params={...})` then
  `submit_job`. This is where the GPU job runs.

They share the workspace directory but not memory — hand the build result over
through a JSON file.

## Build once

```python
# compute_provider tool, provider='modal'
r = build_env('proteomics_gpu', hydrate=True)     # or 'proteomics_boltz_gpu'
import json; json.dump(r, open('built_env.json', 'w'))
print(r['image'], r['spec_sha'], r['volumes'])
```

Record the result as a `### env:<name>@<spec_sha>` block in `compute_details` so
the next session reads the `im-…` straight from the ledger and skips this step.

## Submit the prediction

Stage the FASTA with `select_structure_model.py` first — it fixes the chain
order that `scan_surface_patches.py` depends on.

```python
# repl tool — cell 1 submits and returns
import json
r = json.load(open('built_env.json'))
c = host.compute.create('modal', provider_params={
    'image':   r['image'],          # the im-… string, not the env name
    'env':     'proteomics_gpu',    # NAME -> the shipped egress declaration applies
    'gpu':     'A100-80GB',
    'volumes': r['volumes'],
})
job = c.submit_job(
    intent='Chai-1 on Ab1 Fv (VH+VL), apo',
    command='mkdir -p out && chai-lab fold Ab1_Fv.fasta out/',
    inputs=['Ab1_Fv.fasta'],
    outputs=['*.cif', '*.json'],
    run_timeout_s=1800,
)
```

Then exit the cell and `wait_for_notification`. Only what the job writes under
`./out/` comes back. Read the record afterwards and close the handle:

```python
# repl tool — cell 2, after compute_done
res = c.attach_job('<JOB_ID>').result()
print(res.stdout_tail)
c.close(intent='prediction done — stop the sandbox')
```

For Boltz, swap the env/image and the command (`boltz predict …`); the
submit shape is identical.

## Back to the report

The harvested models land under `hpc/<job_id>/`. Select the best by aggregate
score (0.2×pTM + 0.8×ipTM), then:

```bash
uv run scripts/scan_surface_patches.py hpc/<job_id>/pred.model_idx_0.cif \
    --chains A B --cdr-json cdrmap_chains.json --json > patches.json
uv run scripts/compile_report.py ... --patches patches.json > report.md
```

That fills **Phase 6 — Structural Developability**, which otherwise reports
blocked. Chain IDs must be the ones `select_structure_model.py` assigned.

## Cost discipline

- One handle is one container; it bills idle until `c.close()`, 15 minutes of
  inactivity, or its container timeout.
- Close every handle that ran a job — except one holding unrecovered outputs
  (`result_rejected` / `harvest_failed`), where the sandbox has the only copy.
- `host.compute.ledger()` is the cheap final check that nothing is still live.
