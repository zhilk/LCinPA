function results = fit_lcm_simple(subT, scoreMask)
%FIT_LCM_SIMPLE  Simple three-cause latent-cause model with perceptual
%  read-out.
%
%      R(t) = sum_k post_full(k) * mu_US(k)    (posterior given CS and US)
%
%  Three fixed slots per block, count prior with empty slots at c0 = 0.1,
%  causes reset at block boundaries. sigma_us (pain informativeness for
%  cause assignment) is the only fitted model parameter, grid maximum
%  likelihood over 12 values; sigma_cs = 0.2 and pseudocount a = 0.1 fixed
%  (sweeps in results/). Learning on all trials; likelihood scored on
%  valid test trials unless scoreMask overrides.
%  nPar = 4 (sigma_us, b0, b1, sigma).

    sub = prep_subject(subT);
    if nargin < 2, scoreMask = sub.vt; end
    g = model_grids();

    [pred, params] = predict_model(sub, 'simple', g);
    f = fit_on_grid(pred, params, sub.CR, scoreMask, g.npar.simple);

    results.sigma_us = f.p(1);
    results.sigma_cs = 0.2;          % fixed
    results.c0       = 0.1;          % fixed
    results.a        = 0.1;          % fixed
    results.mu0      = 0.5;
    results.R        = pred(:, f.idx);
    results.CRpred   = f.b0 + f.b1*results.R;
    results.LL       = f.LL;
    results.nPar     = g.npar.simple;
    results.BIC      = f.BIC;
    results.b0 = f.b0; results.b1 = f.b1; results.sigma = f.sigma;
    results.scoreMask = scoreMask;

    % latent structure at the ML parameter: per-trial assignments with
    % cause ids unique across blocks, causes founded during conditioning,
    % and the fraction of morph trials assigned to a new cause
    z = zeros(sub.n, 1); condCauses = []; offset = 0;
    for bi = 1:numel(sub.ub)
        ix = find(sub.blocks == sub.ub(bi));
        r  = lcm_infer_simple(sub.X(ix, :), results.sigma_us);
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
