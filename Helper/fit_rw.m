function results = fit_rw(subT)
%  Fit Rescorla-Wagner model to one participant's trial-level data.

%  Extension of the BCC (Büchel et al., 2014) model with in which perceived pain is a
%  precision-weighted combination of a  prior and the bottom-up sensory
%  signal:
%
%      CR = kappa * V(t) + (1 - kappa) * s(t)
% 
%  where the expectation V(t) is the cue vector CS with learned
%  weights updated by the Rescorla Wagner delta updating rule. 
%       
%      V(t) = CS(t) * w(t)
%      w(t) = x(t-1) + eta * (s(t-1) - V(t-1)) * CS(t)
%

%
%  INPUT:
%    subT – table for one subject, sorted by TrialGlobal, with columns:
%         x_face, x_house, TargetVAS, VASRating, Phase, Block, VisualCategory
%
%  OUTPUT:
%    results – struct with fields:
%      .kappa     – weight parameter
%      .eta       – best-fit learning rate
%      .V         – [nTrials x 1] model predictions (raw, before rescaling)
%      .CRpred    – [nTrials x 1] rescaled predictions
%      .LL        – maximized log-likelihood
%      .BIC       – Bayesian Information Criterion
%      .nPar      – effective number of free parameters (eta + b0, b1, sigma)
%      .b0, .b1   – rescaling parameters
%      .sigma     – noise SD

    subT = sortrows(subT, 'TrialGlobal');

    % --- per-block conditioning w0_seq ---
    blocks = subT.Block;
    if iscell(blocks); blocks = string(blocks); end
    ub = unique(blocks, 'stable');          % block labels in order

    w0_byblock = cell(numel(ub),1);
    for bi = 1:numel(ub)
        c = subT(blocks==ub(bi) & strcmp(subT.Phase,'conditioning'), :);
        avh = mean(c.VASRating(strcmp(c.VisualCategory,'house')), 'omitnan');
        avf = mean(c.VASRating(strcmp(c.VisualCategory,'face')),  'omitnan');
        w0_byblock{bi} = [avf; avh]/100;
    end

    % test trials, tagged by which block they belong to ---
    testMask = ~strcmp(subT.Phase,'conditioning'); % & subT.VASResponse==1 & subT.CatchTrial==0 & ~isnan(subT.VASRating);
    testT    = subT(testMask,:);
    tb       = testT.Block; if iscell(tb); tb = string(tb); end

    % map w0 to its block's 
    w0_seq = zeros(2, height(testT));       % 2 features x nTestTrials
    for bi = 1:numel(ub)
        w0_seq(:, tb==ub(bi)) = repmat(w0_byblock{bi}, 1, nnz(tb==ub(bi)));
    end

    % reset index = first test trial of each block after the first
    resetIdx = find(tb(2:end) ~= tb(1:end-1)) + 1;

    CS = [testT.x_face, testT.x_house];
    US = testT.TargetVAS / 100;
    CR = testT.VASRating;
    nTrials = height(testT);

    % Grid search over eta
    etaGrid = linspace(0.01, 1, 50);
    kappaGrid = linspace(0, 1, 30);
    bestLL  = -Inf;
    bestEta = NaN;
    bestKappa = NaN;
    bestR   = [];

    for k = 1:numel (kappaGrid)
        for i = 1:numel(etaGrid)
            eta = etaGrid(i);
            kappa = kappaGrid (k);
            [R, ~]   = rw_forward(CS, US, eta, kappa, w0_seq, resetIdx);
            LL  = rescaled_LL(R, CR);
            if LL > bestLL
                bestLL  = LL;
                bestEta = eta;
                bestKappa = kappa; 
                bestR   = R;
            end
        end
    end

    % Refine with fmincon
    obj = @(p) -rescaled_LL(rw_forward(CS, US, p(1), p(2), w0_seq, resetIdx), CR);
    opts_opt = optimoptions('fmincon','Display','off');
    pOpt = fmincon(obj, [bestEta bestKappa], [],[],[],[], [0.001 0], [1 1], [], opts_opt);
    etaOpt   = pOpt(1);kappaOpt = pOpt(2);

    [R, W]  = rw_forward(CS, US, etaOpt, kappaOpt, w0_seq, resetIdx);
    [LL, b0, b1, sigma] = rescaled_LL(R, CR);

    % Store results
    results.w0     = w0_byblock;
    results.w      = W;
    results.eta    = etaOpt;
    results.kappa  = kappaOpt; 
    results.R      = R;
    results.CRpred = b0 + b1 * R;
    results.LL     = LL;
    results.nPar   = 5;                   % eta, kappa, b0, b1, sigma
    results.BIC    = -2*LL + results.nPar * log(nTrials);
    results.b0     = b0;
    results.b1     = b1;
    results.sigma  = sigma;
end


%% ========================================================================
function [R, W] = rw_forward(CS, US, eta, kappa, w0_seq, resetIdx)
    nTrials = size(CS, 1);
    w = w0_seq(:,1);                       % initialize at av rating from conditioning
    V = zeros(nTrials, 1);
    R = zeros(nTrials, 1);
    W = zeros(nTrials, numel(w));          % weight vector 
    for t = 1:nTrials
        if any(t == resetIdx)
            w = w0_seq(:,t);           % reset to this block's conditioning prior
        end
        W(t,:) = w';
        x = CS(t,:)';
        s = US(t); 
        V(t) = w' * x;
        R(t) = kappa * V(t) + (1-kappa) * s;
        % weight updating
        delta = US(t) - V(t);
        w = w + eta * x * delta;
    end
end


%% ========================================================================
function [LL, b0, b1, sigma] = rescaled_LL(R, CR)
    n  = numel(CR);
    X  = [ones(n,1), R(:)];
    b  = X \ CR(:);                       % OLS
    b0 = b(1);
    b1 = b(2);
    pred  = X * b;
    resid = CR(:) - pred;
    sigma = sqrt(mean(resid.^2));         % ML estimate of sigma; assumes homoskedasticity
    LL    = -0.5 * n * log(2*pi) - n*log(sigma) - 0.5*sum((resid/sigma).^2);
end
