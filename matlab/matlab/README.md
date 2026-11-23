# MATLAB publication code: models, fits, parameter recovery, model recovery

Self-contained MATLAB implementation of the paper's model set and the
recovery analyses. Everything here is scored with one convention and the
recovery scripts reuse the same forward passes as the fitters. The older
`models/` folder (coauthor originals and earlier packaged fitters) is kept
unchanged as provenance; this folder is the code of record.

## Model set

| fitter | model | free params | grid | nPar |
|---|---|---|---|---|
| `fit_null.m` | R = US (no expectation) | - | - | 3 |
| `fit_rw.m` | Rescorla-Wagner: R = w'CS | eta | 50 (0.01-1) | 4 |
| `fit_bcc.m` | Bayesian cue combination, fixed objective weights: R = kappa*V + (1-kappa)*US, V = CS * w_obj | kappa | 30 (0-1) | 4 |
| `fit_rwbcc.m` | RW learning + Bayesian integration: R = kappa*(w'CS) + (1-kappa)*US | eta, kappa | 30 x 30 | 5 |
| `fit_lcm_simple.m` | 3-cause LCM, count prior (c0 = 0.1) | sigma_us | 12 values | 4 |
| `fit_lcm_percept.m` | CRP LCM, perceptual read-out | alpha, sigma_us | 12 x 12 | 5 |

nPar counts include the rating-model parameters b0, b1, sigma. The full
LCM (7 parameters) is excluded from the model set; its identifiability
audit is reported in the supplement.

## Conventions (identical across all six models)


- **Stimulus coding.** Stimuli are coded as [1 0] house, [0 1] face and [0.5 0.5] for the morph.
- **Per-block reset.** Every model begins each block naive: RW weights
  return to [0.5 0.5], the BCC objective weights are defined per block, and
  latent causes form fresh within each block. Motivated by the design: the
  cue-outcome contingencies are re-conditioned in every block.
- **Learning on all trials, likelihood scored on valid test trials**
  (rating given, conditioning phase excluded). Fitters accept an optional
  second argument to score a different trial mask.
- **Rating model.** Each model produces a per-trial internal prediction
  R(t) from the stimuli alone; ratings never enter the forward pass.
  Ratings are scored as CR = b0 + b1*R + N(0, sigma^2), solved in closed
  form (`rescaled_LL.m` / vectorized in `fit_on_grid.m`).
- **Grid maximum likelihood** (no gradient refinement), grids above.
- **BIC** = -2*LL + nPar*log(n), n = number of scored trials.

## Files

Fitters (`fit_*.m`): one subject table in, one results struct out
(parameters, LL, BIC, b0, b1, sigma, R, CRpred; the LCM fitters add the
latent assignment structure).

Shared machinery:
- `load_subjects.m` / `prep_subject.m` - data loading (studies 1-3,
  N = 156; study 3 subject 36 excluded) and per-subject trial masks
- `model_grids.m` - all parameter grids and nPar in one place
- `predict_model.m` - per-subject prediction matrix, one column per grid
  point (predictions depend on stimuli only, so they are computed once and
  reused across real fits and every synthetic refit)
- `fit_on_grid.m` - vectorized closed-form rating fit over all grid columns
- `rw_forward.m`, `bcc_fixed_V.m`, `lcm_infer_simple.m`,
  `lcm_infer_percept.m` - the forward passes
- `rescaled_LL.m` - reference scalar implementation of the rating model
- `eff_pred.m`, `phenotype_simple.m` - derived read-outs (predicted morph
  effect; learned-vs-novel assignment phenotype)

Analyses:
- `fit_all_subjects.m` - fits all six models to every subject;
  writes `results/matlab_fits.csv`
- `parameter_recovery.m` - 10 synthetic datasets per subject per model,
  generated at the subject's own ML parameters on the subject's own trial
  sequence, refit with the identical pipeline; writes
  `results/recovery_param_simple.csv`, `recovery_param_percept.csv`,
  `recovery_param_rwbcc.csv` (rw, bcc, rwbcc). Seed rng(7).
- `model_recovery.m` - 3 synthetic datasets per subject per generating
  model; all six models refit to each; records the BIC winner (valid test
  trials) and the morph-trial LL winner (at the valid-test ML parameters,
  no refit); writes `results/recovery_model_confusion.csv`. Seed rng(11).

## Running

With this folder on the path (data paths resolve from the folder's own
location, so the working directory does not matter):

```matlab
fit_all_subjects
parameter_recovery
model_recovery
```

Scripts use plain loops (no parfor) so the seeded runs are exactly
reproducible. Do not put `models/` and `matlab/` on the path together; the
fitter names overlap.

## Robustness note

All code and all reported results use the per-block reset. As a
sensitivity check, RW and RW-BCC were additionally fit with weights
carried across blocks (the convention of some earlier analyses); the
alternative changes little: RW median per-subject delta LL = 0.01,
summed -11; RW-BCC summed +110 under carry, with kappa medians 0.45
(carry) vs 0.26 (reset); 13/156 BIC winners flip; no conclusion changes.
