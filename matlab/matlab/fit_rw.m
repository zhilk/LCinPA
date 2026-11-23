function results = fit_rw(subT, scoreMask)
%FIT_RW  Rescorla-Wagner model: ratings driven by the learned cue
%  expectation alone.
%
%      R(t)    = w' * CS(t)
%      w_{t+1} = w + eta * CS(t) * (US(t) - R(t))
%
%  Weights start naive at [0.5 0.5] and reset at each block boundary. 
%  Eta fit by grid maximum likelihood (50 values, 0.01-1).
%  Note here that R(t) = V(t) as the rating prediction is equivalent to the
%  expectation for this model. 

%  Learning on all trials; likelihood scored on valid test trials unless
%  scoreMask overrides. 
% 
%  nPar = 4 (eta, b0, b1, sigma).

    sub = prep_subject(subT);
    if nargin < 2, scoreMask = sub.vt; end
    g = model_grids();

    [pred, params] = predict_model(sub, 'rw', g);
    f = fit_on_grid(pred, params, sub.CR, scoreMask, g.npar.rw);

    results.eta    = f.p(1);
    results.R      = pred(:, f.idx);
    results.CRpred = f.b0 + f.b1*results.R;
    results.LL     = f.LL;
    results.nPar   = g.npar.rw;
    results.BIC    = f.BIC;
    results.b0 = f.b0; results.b1 = f.b1; results.sigma = f.sigma;
    results.w_0 = [0.5; 0.5];
    results.scoreMask = scoreMask;
end
