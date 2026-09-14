function results = fit_rw(T)
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
%    T – table for one subject, sorted by TrialGlobal, with columns:
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

    T = sortrows(T, 'TrialGlobal');

    % determine w_fixed 
    keep = strcmp(T.Phase,'conditioning'); %& T.VASResponse==1 & T.CatchTrial==0;
    C = T(keep,:);
    av_house = mean(C.VASRating(strcmp(C.VisualCategory,'house')), 'omitnan');
    av_face  = mean(C.VASRating(strcmp(C.VisualCategory,'face')),  'omitnan');
    w_0 = [av_face; av_house]/100;

    % extract experimental data from test phase 
    testMask = ~strcmp(T.Phase,'conditioning'); % & T.VASResponse==1 & T.CatchTrial==0 & ~isnan(T.VASRating);
    testT = T(testMask,:);

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
            R   = rw_forward(CS, US, eta, kappa, w_0);
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
    obj = @(p) -rescaled_LL(rw_forward(CS, US, p(1), p(2), w_0), CR);
    opts_opt = optimoptions('fmincon','Display','off');
    pOpt = fmincon(obj, [bestEta bestKappa], [],[],[],[], [0.001 0], [1 1], [], opts_opt);
    etaOpt   = pOpt(1);kappaOpt = pOpt(2);

    R  = rw_forward(CS, US, etaOpt, kappaOpt, w_0);
    [LL, b0, b1, sigma] = rescaled_LL(R, CR);

    % Store results
    results.w0     = w0;
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
function R = rw_forward(CS, US, eta, kappa, w_0)
    nTrials = size(CS, 1);
    w = w_0;                       % initialize at midpoint
    V = zeros(nTrials, 1);
    R = zeros(nTrials, 1);
    for t = 1:nTrials
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
