function results = fit_rwbcc(subT, scoreMask)
%FIT_RWBCC  RW learning combined with Bayesian integration: the learned
%  cue expectation is blended with the delivered stimulus.
%
%      V(t)   = w' * CS(t)                            (RW delta rule, see fit_rw for weight updating)
%      w_{t+1} = w + eta * CS(t) * (US(t) - V(t))     (rw forward)
%      R(t)   = kappa * V(t) + (1 - kappa) * US(t)
%
%  kappa = 1 reduces to the pure RW model, kappa = 0 to the null model.
%  Weights reset naive at each block boundary. 
% 
%  eta and kappa fit jointly by grid maximum likelihood (30 x 30). 
%  Learning on all trials;likelihood scored on valid test trials unless scoreMask overrides.
%  
%  nPar = 5 (eta, kappa, b0, b1, sigma).

    sub = prep_subject(subT);
    if nargin < 2, scoreMask = sub.vt; end
    g = model_grids();

    [pred, params] = predict_model(sub, 'rwbcc', g);
    f = fit_on_grid(pred, params, sub.CR, scoreMask, g.npar.rwbcc);

    results.eta    = f.p(1);
    results.kappa  = f.p(2);
    results.R      = pred(:, f.idx);
    results.CRpred = f.b0 + f.b1*results.R;
    results.LL     = f.LL;
    results.nPar   = g.npar.rwbcc;
    results.BIC    = f.BIC;
    results.b0 = f.b0; results.b1 = f.b1; results.sigma = f.sigma;
    results.w_0 = [0.5; 0.5];
    results.scoreMask = scoreMask;
end
