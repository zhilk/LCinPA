function [pred, params] = predict_model(sub, model, g)
%PREDICT_MODEL  Prediction matrix of one model for one subject: one column
%  per grid point.
%
%  All model predictions R(t) depend on the stimuli only (never on the
%  ratings), so each grid point's forward pass is computed once per subject
%  and reused for the real fit and for every synthetic refit in the
%  recovery analyses - fitting a dataset then reduces to the closed-form
%  rating fit per column (fit_on_grid.m).
%
%  pred    [nTrials x nGrid] prediction matrix
%  params  [nGrid x nFree] model parameter value(s) per column:
%    'null'     []                    (single column, R = US)
%    'rw'       eta                   (50)
%    'bcc'      kappa                 (30)       R = k*V_obj + (1-k)*US
%    'rwbcc'    [eta kappa]           (30 x 30)  R = k*mu_RW + (1-k)*US
%    'simple'   sigma_us              (12)
%    'percept'  [alpha sigma_us]      (12 x 12)
%  Column order fixes the tie-breaking convention (first column of equal
%  likelihood wins); for 'percept', alpha varies slowest. At alpha = 0 the
%  prediction is sigma_us-independent, so ties resolve to the smallest
%  sigma_us by this ordering.

    switch model
        case 'null'
            pred = sub.US;
            params = zeros(1, 0);

        case 'rw'
            nG = numel(g.eta);
            pred = zeros(sub.n, nG);
            for i = 1:nG
                pred(:, i) = rw_forward(sub.CS, sub.US, g.eta(i), sub.blocks);
            end
            params = g.eta(:);

        case 'bcc'
            V = bcc_fixed_V(sub);
            nG = numel(g.kappa);
            pred = zeros(sub.n, nG);
            for i = 1:nG
                pred(:, i) = g.kappa(i)*V + (1 - g.kappa(i))*sub.US;
            end
            params = g.kappa(:);

        case 'rwbcc'
            nE = numel(g.eta30); nK = numel(g.kappa);
            mus = zeros(sub.n, nE);
            for i = 1:nE
                mus(:, i) = rw_forward(sub.CS, sub.US, g.eta30(i), sub.blocks);
            end
            pred = zeros(sub.n, nE*nK);
            params = zeros(nE*nK, 2);
            c = 0;
            for i = 1:nE
                for j = 1:nK
                    c = c + 1;
                    pred(:, c) = g.kappa(j)*mus(:, i) + (1 - g.kappa(j))*sub.US;
                    params(c, :) = [g.eta30(i), g.kappa(j)];
                end
            end

        case 'simple'
            nG = numel(g.sus);
            pred = zeros(sub.n, nG);
            for i = 1:nG
                for bi = 1:numel(sub.ub)
                    m = sub.blocks == sub.ub(bi);
                    r = lcm_infer_simple(sub.X(m, :), g.sus(i));
                    pred(m, i) = r.R;
                end
            end
            params = g.sus(:);

        case 'percept'
            nA = numel(g.alpha); nS = numel(g.sus);
            pred = zeros(sub.n, nA*nS);
            params = zeros(nA*nS, 2);
            c = 0;
            for ia = 1:nA
                for is = 1:nS
                    c = c + 1;
                    opts = struct('alpha', g.alpha(ia), 'stickiness', 0, ...
                                  'a', 0.1, 'mu0', 0.5, 'K', 10, ...
                                  'sigma_cs', 0.2, 'sigma_us', g.sus(is));
                    for bi = 1:numel(sub.ub)
                        m = sub.blocks == sub.ub(bi);
                        r = lcm_infer_percept(sub.X(m, :), opts);
                        pred(m, c) = r.R;
                    end
                    params(c, :) = [g.alpha(ia), g.sus(is)];
                end
            end

        otherwise
            error('unknown model %s', model);
    end
end
