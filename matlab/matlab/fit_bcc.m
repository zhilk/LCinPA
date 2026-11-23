function results = fit_bcc(subT, scoreMask)
%FIT_BCC  Bayesian cue combination with fixed objective weights: no
%  learning, the cue expectation is normatively correct from trial one.
%
%      V(t) = CS(t) * w_obj(block)      (w_obj = delivered heat per cue during conditioning)
%      R(t) = kappa * V(t) + (1 - kappa) * US(t)
%
%  kappa arbitrates between the cue-based expectation and the delivered
%  stimulus; kappa = 0 reduces to the null model. 
% 
%  kappa fit by grid maximum likelihood (30 values, 0-1). 
%  
%  Likelihood scored on valid test trials unless scoreMask overrides. 
% 
%  nPar = 4 (kappa, b0, b1, sigma).

    sub = prep_subject(subT);
    if nargin < 2, scoreMask = sub.vt; end
    g = model_grids();

    [pred, params] = predict_model(sub, 'bcc', g);
    f = fit_on_grid(pred, params, sub.CR, scoreMask, g.npar.bcc);

    results.kappa  = f.p(1);
    results.V      = bcc_fixed_V(sub);
    results.R      = pred(:, f.idx);
    results.CRpred = f.b0 + f.b1*results.R;
    results.LL     = f.LL;
    results.nPar   = g.npar.bcc;
    results.BIC    = f.BIC;
    results.b0 = f.b0; results.b1 = f.b1; results.sigma = f.sigma;
    results.scoreMask = scoreMask;
end
