function results = fit_lcm_percept(subT, scoreMask)
%FIT_LCM_PERCEPT  CRP latent-cause model with perceptual read-out.
%
%      R(t) = sum_k post_full(k) * mu_US(k)    (posterior given CS and US)
%
%  Causes form under a CRP prior (concentration alpha, capacity K = 10)
%  within each block and reset at block boundaries. alpha and sigma_us fit
%  jointly by grid maximum likelihood (12 x 12; the alpha grid reaches 50
%  because values above ~10 are behaviorally saturated and fits pile up at
%  a lower ceiling otherwise). alpha = 0 reduces to a single running
%  cause. sigma_cs = 0.2 and pseudocount a = 0.1 fixed (see
%  lcm_infer_percept.m for the role of a). Learning on all trials;
%  likelihood scored on valid test trials unless scoreMask overrides.
%  nPar = 5 (alpha, sigma_us, b0, b1, sigma).

    sub = prep_subject(subT);
    if nargin < 2, scoreMask = sub.vt; end
    g = model_grids();

    [pred, params] = predict_model(sub, 'percept', g);
    f = fit_on_grid(pred, params, sub.CR, scoreMask, g.npar.percept);

    results.alpha    = f.p(1);
    results.sigma_us = f.p(2);
    results.sigma_cs = 0.2;          % fixed
    results.a        = 0.1;          % fixed
    results.mu0      = 0.5;
    results.R        = pred(:, f.idx);
    results.CRpred   = f.b0 + f.b1*results.R;
    results.LL       = f.LL;
    results.nPar     = g.npar.percept;
    results.BIC      = f.BIC;
    results.b0 = f.b0; results.b1 = f.b1; results.sigma = f.sigma;
    results.scoreMask = scoreMask;

    % alpha = 0 baseline: best single-cause fit over the sigma_us grid
    % (identical columns; the first alpha = 0 column represents them)
    i0 = find(params(:,1) == 0, 1);
    f0 = fit_on_grid(pred(:, i0), params(i0, :), sub.CR, scoreMask, g.npar.percept);
    results.LL_single = f0.LL;
    results.logLR_vs_single = f.LL - f0.LL;

    % latent structure at the ML parameters (as in fit_lcm_simple)
    opts = struct('alpha', results.alpha, 'stickiness', 0, 'a', 0.1, ...
                  'mu0', 0.5, 'K', 10, 'sigma_cs', 0.2, ...
                  'sigma_us', results.sigma_us);
    z = zeros(sub.n, 1); condCauses = []; offset = 0;
    for bi = 1:numel(sub.ub)
        ix = find(sub.blocks == sub.ub(bi));
        r  = lcm_infer_percept(sub.X(ix, :), opts);
        act   = find(r.Nk > 0);
        remap = zeros(1, numel(r.Nk));
        remap(act) = 1:numel(act);
        z(ix) = remap(r.z)' + offset;
        ci = sub.isCond(ix);
        for j = 1:numel(act)
            fk = find(r.z == act(j), 1);          % first trial in this cause
            if ~isempty(fk) && ci(fk)
                condCauses(end+1) = j + offset; %#ok<AGROW>
            end
        end
        offset = offset + numel(act);
    end
    results.z = z;
    results.condCauses = condCauses;
    results.morph_z = z(sub.isMix);
    results.morph_newcause_frac = mean(~ismember(z(sub.isMix), condCauses));
end
