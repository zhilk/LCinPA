function results = fit_null(subT, scoreMask)
%FIT_NULL  Null model: ratings driven by the delivered stimulus alone.
%
%      R(t) = US(t)
%
%  No expectation, no learning, no latent causes; the affine rating model
%  (b0, b1, sigma) is the only fitted part, so this is the baseline any
%  expectation term must beat. Likelihood scored on valid test trials
%  unless scoreMask overrides. nPar = 3 (b0, b1, sigma).

    sub = prep_subject(subT);
    if nargin < 2, scoreMask = sub.vt; end
    g = model_grids();

    [pred, params] = predict_model(sub, 'null', g);
    f = fit_on_grid(pred, params, sub.CR, scoreMask, g.npar.null);

    results.R      = pred;
    results.CRpred = f.b0 + f.b1*pred;
    results.LL     = f.LL;
    results.nPar   = g.npar.null;
    results.BIC    = f.BIC;
    results.b0 = f.b0; results.b1 = f.b1; results.sigma = f.sigma;
    results.scoreMask = scoreMask;
end
