function results = fit_rw(T)
%FIT_RW  Fit Rescorla-Wagner model to one participant's trial-level data.
%
%  INPUT:
%    T – table for one subject, sorted by TrialGlobal, with columns:
%         x_face, x_house, TargetVAS, VASRating, Phase, Block, VisualCategory
%
%  OUTPUT:
%    results – struct with fields:
%      .eta       – best-fit learning rate
%      .V         – [nTrials x 1] model predictions (raw, before rescaling)
%      .CRpred    – [nTrials x 1] rescaled predictions
%      .LL        – maximized log-likelihood
%      .BIC       – Bayesian Information Criterion
%      .nPar      – effective number of free parameters (eta + b0, b1, sigma)
%      .b0, .b1   – rescaling parameters
%      .sigma     – noise SD

    % Sort by global trial order
    T = sortrows(T, 'TrialGlobal');
    nTrials = height(T);

    % Extract CS and US
    CS = [T.x_face, T.x_house];           % [nTrials x 2]
    US = T.TargetVAS / 100;               % normalize to [0,1]
    CR = T.VASRating;                      % raw VAS ratings

    % Grid search over eta
    etaGrid = linspace(0.01, 1, 50);
    bestLL  = -Inf;
    bestEta = NaN;
    bestV   = [];

    for i = 1:numel(etaGrid)
        eta = etaGrid(i);
        V   = rw_forward(CS, US, eta);
        LL  = rescaled_LL(V, CR);
        if LL > bestLL
            bestLL  = LL;
            bestEta = eta;
            bestV   = V;
        end
    end

    % Refine with fmincon
    obj = @(p) -rescaled_LL(rw_forward(CS, US, p), CR);
    opts_opt = optimoptions('fmincon','Display','off');
    etaOpt = fmincon(obj, bestEta, [],[],[],[], 0.001, 1, [], opts_opt); %fminbnd might work here, as it is only one variable

    V  = rw_forward(CS, US, etaOpt);
    [LL, b0, b1, sigma] = rescaled_LL(V, CR);

    % Store results
    results.eta    = etaOpt;
    results.V      = V;
    results.CRpred = b0 + b1 * V;
    results.LL     = LL;
    results.nPar   = 4;                   % eta, b0, b1, sigma
    results.BIC    = -2*LL + results.nPar * log(nTrials);
    results.b0     = b0;
    results.b1     = b1;
    results.sigma  = sigma;
end


%% ========================================================================
function V = rw_forward(CS, US, eta)
%RW_FORWARD  Run RW model forward, return prediction on each trial.
    nTrials = size(CS, 1);
    w = [0.5; 0.5];                       % initialize at midpoint
    V = zeros(nTrials, 1);
    for t = 1:nTrials
        x = CS(t,:)';
        V(t) = w' * x;
        delta = US(t) - V(t);
        w = w + eta * x * delta;
    end
end


%% ========================================================================
function [LL, b0, b1, sigma] = rescaled_LL(V, CR)
%RESCALED_LL  Compute log-likelihood after linear rescaling V -> CR.
    n  = numel(CR);
    X  = [ones(n,1), V(:)];
    b  = X \ CR(:);                       % OLS
    b0 = b(1);
    b1 = b(2);
    pred  = X * b;
    resid = CR(:) - pred;
    sigma = sqrt(mean(resid.^2));         % ML estimate of sigma; assumes homoskedasticity
    LL    = -0.5 * n * log(2*pi) - n*log(sigma) - 0.5*sum((resid/sigma).^2);
end
