function results = fit_lcm(T, resetBlocks)

%   Fit Latent Cause Model to one participant's data.
%
%  Marginalizes over alpha (concentration parameter) via grid integration.
%  Variance parameters are optimized via outer grid search.
%
%  INPUTS:
%    T            – table for one subject (sorted by TrialGlobal)
%    resetBlocks  – logical; if true, reinitialize causes at block boundary
%                   (default: false = full carry-over)
%
%  OUTPUT:
%    results – struct with fields:
%      .alpha_map   – posterior mean of alpha
%      .sigma_cs    – best-fit CS variance
%      .sigma_us    – best-fit US variance
%      .V           – [nTrials x 1] model predictions (raw)
%      .CRpred      – [nTrials x 1] rescaled predictions
%      .LL          – marginal log-likelihood (integrated over alpha)
%      .LL_alpha0   – log-likelihood at alpha=0 (single-cause = nested RW)
%      .logBF_vs_rw – log Bayes Factor: LCM vs single-cause (alpha=0)
%      .BIC         – BIC (using 3 effective params + rescaling)
%      .nPar        – effective parameter count
%      .z           – [nTrials x 1] MAP cause assignments

    if nargin < 2; resetBlocks = false; end

    T = sortrows(T, 'TrialGlobal');
    nTrials = height(T);

    % Build stimulus matrix: [US, CS_face, CS_house]
    US = T.TargetVAS / 100;
    CS = [T.x_face, T.x_house];
    CR = T.VASRating;
    X  = [US, CS];

    % Find block boundary (if two blocks exist)
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

    % --- Grid over variance parameters ---
    sigmaGrid = [0.1 0.25 0.5 0.75 1.0];
    nAlpha    = 50;
    alphaGrid = linspace(0, 10, nAlpha);

    bestMargLL = -Inf;
    bestSigCS  = NaN;
    bestSigUS  = NaN;
    bestAlphaPost = [];
    bestVall   = [];
    bestZall   = [];

    for sc = 1:numel(sigmaGrid)
        for su = 1:numel(sigmaGrid)
            sigma_cs = sigmaGrid(sc);
            sigma_us = sigmaGrid(su);

            % Evaluate LL for each alpha
            logLiks = zeros(1, nAlpha);
            Vall    = cell(1, nAlpha);
            Zall    = cell(1, nAlpha);

            for ia = 1:nAlpha
                alpha = alphaGrid(ia);
                opts  = struct('alpha', alpha, 'stickiness', 0, ...
                               'sigma_cs', sigma_cs, 'sigma_us', sigma_us, ...
                               'a', 1, 'b', 1, 'K', 10);

                if isempty(blockChange)
                    % Full carry-over
                    res = lcm_infer(X, opts);
                else
                    % Reset at block boundary
                    res1 = lcm_infer(X(1:blockChange-1, :), opts);
                    res2 = lcm_infer(X(blockChange:end, :), opts);
                    res.V = [res1.V; res2.V];
                    res.z = [res1.z; res2.z];
                end

                Vall{ia} = res.V;
                Zall{ia} = res.z;
                logLiks(ia) = rescaled_LL(res.V, CR);
            end

            % Marginal log-likelihood: log(1/N * sum(exp(LL)))
            margLL = logsumexp(logLiks) - log(nAlpha);

            if margLL > bestMargLL
                bestMargLL    = margLL;
                bestSigCS     = sigma_cs;
                bestSigUS     = sigma_us;
                bestAlphaPost = logLiks;
                bestVall      = Vall;
                bestZall      = Zall;
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

    % Re-run lcm_infer at MAP alpha to get full result with sufficient stats
    opts_map = struct('alpha', alphaGrid(iMAP), 'stickiness', 0, ...
                      'sigma_cs', bestSigCS, 'sigma_us', bestSigUS, ...
                      'a', 1, 'b', 1, 'K', 10);
    if isempty(blockChange)
        res_map = lcm_infer(X, opts_map);
        res_map.resetUsed = false;
    else
        res_map1 = lcm_infer(X(1:blockChange-1, :), opts_map);
        res_map2 = lcm_infer(X(blockChange:end, :), opts_map);

        % Offset block-2 cause IDs so they don't collide with block-1
        nK1 = numel(res_map1.Nk);
        z2_offset = res_map2.z + nK1;

        % Merge into single result struct
        res_map.V    = [res_map1.V; res_map2.V];
        res_map.z    = [res_map1.z; z2_offset];
        res_map.post = blkdiag(res_map1.post, res_map2.post);
        res_map.K    = res_map1.K + res_map2.K;

        % Concatenate sufficient stats (block-1 causes then block-2 causes)
        res_map.Nk   = [res_map1.Nk, res_map2.Nk];
        res_map.SumF = [res_map1.SumF; res_map2.SumF];
        res_map.mu0  = res_map1.mu0;
        res_map.a_pseudo = res_map1.a_pseudo;

        % Store per-block info for cause-assignment analysis
        res_map.resetUsed   = true;
        res_map.nK_block1   = nK1;
        res_map.nK_block2   = numel(res_map2.Nk);
        res_map.blockChange = blockChange;
    end

    % Posterior-weighted prediction (average over alpha grid)
    V_avg = zeros(nTrials, 1);
    for ia = 1:nAlpha
        V_avg = V_avg + postAlpha(ia) * bestVall{ia};
    end

    % Rescale to VAS
    [LL, b0, b1, sigma] = rescaled_LL(V_avg, CR);

    % LL at alpha = 0 (single-cause, nested within LCM)
    LL_alpha0 = bestAlphaPost(1);  % alphaGrid(1) = 0

    % Store
    results.alpha_map   = alphaMean;
    results.sigma_cs    = bestSigCS;
    results.sigma_us    = bestSigUS;
    results.V           = V_avg;
    results.CRpred      = b0 + b1 * V_avg;
    results.LL          = bestMargLL;
    results.LL_alpha0   = LL_alpha0;
    results.logBF_vs_rw = bestMargLL - LL_alpha0;
    results.nPar        = 6;              % sigma_cs, sigma_us, alpha (marginal) + b0, b1, sigma
    results.BIC         = -2*bestMargLL + results.nPar * log(nTrials);
    results.b0          = b0;
    results.b1          = b1;
    results.sigma       = sigma;
    results.z           = z_best;
    results.alphaGrid   = alphaGrid;
    results.alphaPost   = postAlpha;
    results.Nk          = res_map.Nk;     % cause counts at end of task
    results.SumF        = res_map.SumF;   % sum of features per cause
    results.mu0         = res_map.mu0;
    results.a_pseudo    = res_map.a_pseudo;
    results.post_full   = res_map.post;   % full posterior per trial
    results.resetUsed   = res_map.resetUsed;
    if res_map.resetUsed
        results.nK_block1   = res_map.nK_block1;
        results.nK_block2   = res_map.nK_block2;
        results.blockChange = res_map.blockChange;
    end
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
