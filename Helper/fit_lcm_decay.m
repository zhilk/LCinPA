function results = fit_lcm_decay(T, resetBlocks)
%FIT_LCM_DECAY  Fit LCM with temporal decay to one participant's data.
%
%  Same as fit_lcm but adds a grid search over the decay parameter lambda.
%  Nk = exp(-lambda) * Nk on each trial (Heesink et al., 2024).
%  Lambda = 0 reduces to the standard LCM.
%
%  INPUTS:
%    T            – table for one subject (sorted by TrialGlobal)
%    resetBlocks  – logical; if true, reinitialize causes at block boundary
%                   (default: false = full carry-over)
%
%  OUTPUT:
%    results – struct (same fields as fit_lcm, plus .lambda)

    if nargin < 2; resetBlocks = false; end

    T = sortrows(T, 'TrialGlobal');
    nTrials = height(T);

    % Build stimulus matrix: [US, CS_face, CS_house]
    US = T.TargetVAS / 100;
    CS = [T.x_face, T.x_house];
    CR = T.VASRating;
    X  = [US, CS];

    % Find block boundary
    blocks = T.Block;
    if iscell(blocks); blocks = string(blocks); end
    blockChange = [];
    if resetBlocks
        for t = 2:nTrials
            if ~strcmp(blocks(t), blocks(t-1))
                blockChange = t;
                break
            end
        end
    end

    % --- Parameter grids ---
    sigmaGrid  = [0.1 0.25 0.5 0.75 1.0];
    lambdaGrid = [0 0.01 0.02 0.05 0.1 0.2 0.5];
    nAlpha     = 50;
    alphaGrid  = linspace(0, 10, nAlpha);

    bestMargLL    = -Inf;
    bestSigCS     = NaN;
    bestSigUS     = NaN;
    bestLambda    = NaN;
    bestAlphaPost = [];
    bestVall      = [];
    bestZall      = [];

    for sl = 1:numel(lambdaGrid)
        lam = lambdaGrid(sl);

        for sc = 1:numel(sigmaGrid)
            for su = 1:numel(sigmaGrid)
                sigma_cs = sigmaGrid(sc);
                sigma_us = sigmaGrid(su);

                logLiks = zeros(1, nAlpha);
                Vall    = cell(1, nAlpha);
                Zall    = cell(1, nAlpha);

                for ia = 1:nAlpha
                    alpha = alphaGrid(ia);
                    opts  = struct('alpha', alpha, 'stickiness', 0, ...
                                   'lambda', lam, ...
                                   'sigma_cs', sigma_cs, 'sigma_us', sigma_us, ...
                                   'a', 1, 'b', 1, 'K', 10);

                    if isempty(blockChange)
                        res = lcm_infer(X, opts);
                    else
                        res1 = lcm_infer(X(1:blockChange-1, :), opts);
                        res2 = lcm_infer(X(blockChange:end, :), opts);
                        res.V = [res1.V; res2.V];
                        res.z = [res1.z; res2.z];
                    end

                    Vall{ia} = res.V;
                    Zall{ia} = res.z;
                    logLiks(ia) = rescaled_LL(res.V, CR);
                end

                margLL = logsumexp(logLiks) - log(nAlpha);

                if margLL > bestMargLL
                    bestMargLL    = margLL;
                    bestSigCS     = sigma_cs;
                    bestSigUS     = sigma_us;
                    bestLambda    = lam;
                    bestAlphaPost = logLiks;
                    bestVall      = Vall;
                    bestZall      = Zall;
                end
            end
        end
    end

    % Posterior over alpha
    logPost = bestAlphaPost(:)' - logsumexp(bestAlphaPost);
    postAlpha = exp(logPost);
    alphaMean = sum(postAlpha(:)' .* alphaGrid(:)');

    % Best single alpha (MAP)
    [~, iMAP] = max(bestAlphaPost);
    V_best = bestVall{iMAP};
    z_best = bestZall{iMAP};

    % Posterior-weighted prediction
    V_avg = zeros(nTrials, 1);
    for ia = 1:nAlpha
        V_avg = V_avg + postAlpha(ia) * bestVall{ia};
    end

    [LL, b0, b1, sigma] = rescaled_LL(V_avg, CR);

    % LL at alpha = 0
    LL_alpha0 = bestAlphaPost(1);

    % Store
    results.alpha_map   = alphaMean;
    results.lambda      = bestLambda;
    results.sigma_cs    = bestSigCS;
    results.sigma_us    = bestSigUS;
    results.V           = V_avg;
    results.CRpred      = b0 + b1 * V_avg;
    results.LL          = bestMargLL;
    results.LL_alpha0   = LL_alpha0;
    results.logBF_vs_rw = bestMargLL - LL_alpha0;
    results.nPar        = 7;    % lambda, sigma_cs, sigma_us, alpha (marginal) + b0, b1, sigma
    results.BIC         = -2*bestMargLL + results.nPar * log(nTrials);
    results.b0          = b0;
    results.b1          = b1;
    results.sigma       = sigma;
    results.z           = z_best;
    results.alphaGrid   = alphaGrid;
    results.alphaPost   = postAlpha;
end


%% ========================================================================
function [LL, b0, b1, sigma] = rescaled_LL(V, CR)
    n  = numel(CR);
    X  = [ones(n,1), V(:)];
    b  = X \ CR(:);
    b0 = b(1);  b1 = b(2);
    pred  = X * b;
    resid = CR(:) - pred;
    sigma = sqrt(mean(resid.^2));
    LL    = -0.5*n*log(2*pi) - n*log(sigma) - 0.5*sum((resid/sigma).^2);
end


%% ========================================================================
function s = logsumexp(x)
    mx = max(x);
    s  = mx + log(sum(exp(x - mx)));
end
