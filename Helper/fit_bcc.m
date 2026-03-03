function results = fit_bcc(T)
%FIT_BCC  Fit Bayesian Cue Combination model to one participant's data.
%
%  The BCC model (Büchel et al., 2014) computes perceived pain as a
%  precision-weighted combination of a learned prior (via delta-rule) and
%  the bottom-up sensory signal:
%     CR = kappa * mu_prior + (1 - kappa) * s
%
%  INPUT:
%    T – table for one subject, sorted by TrialGlobal, with columns:
%         x_face, x_house, TargetVAS, VASRating, Phase, Block, VisualCategory
%
%  OUTPUT:
%    results – struct with fields:
%      .eta, .kappa  – best-fit parameters
%      .V            – [nTrials x 1] model predictions (raw)
%      .CRpred       – [nTrials x 1] rescaled predictions
%      .LL, .BIC     – fit statistics
%      .nPar         – effective free parameters (eta, kappa, b0, b1, sigma)

    T = sortrows(T, 'TrialGlobal');
    nTrials = height(T);

    CS = [T.x_face, T.x_house];
    US = T.TargetVAS / 100;
    CR = T.VASRating;

    % Grid search over eta and kappa
    etaGrid   = linspace(0.01, 1, 30);
    kappaGrid = linspace(0, 1, 30);
    bestLL    = -Inf;
    bestP     = [NaN NaN];

    for i = 1:numel(etaGrid)
        for j = 1:numel(kappaGrid)
            V  = bcc_forward(CS, US, etaGrid(i), kappaGrid(j));
            LL = rescaled_LL(V, CR);
            if LL > bestLL
                bestLL = LL;
                bestP  = [etaGrid(i), kappaGrid(j)];
            end
        end
    end

    % Refine with fmincon
    obj = @(p) -rescaled_LL(bcc_forward(CS, US, p(1), p(2)), CR);
    opts_opt = optimoptions('fmincon','Display','off');
    pOpt = fmincon(obj, bestP, [],[],[],[], [0.001 0], [1 1], [], opts_opt);

    V = bcc_forward(CS, US, pOpt(1), pOpt(2));
    [LL, b0, b1, sigma] = rescaled_LL(V, CR);

    results.eta    = pOpt(1);
    results.kappa  = pOpt(2);
    results.V      = V;
    results.CRpred = b0 + b1 * V;
    results.LL     = LL;
    results.nPar   = 5;                    % eta, kappa, b0, b1, sigma
    results.BIC    = -2*LL + results.nPar * log(nTrials);
    results.b0     = b0;
    results.b1     = b1;
    results.sigma  = sigma;
end


%% ========================================================================
function V = bcc_forward(CS, US, eta, kappa)
%BCC_FORWARD  Run BCC model forward.
%  Prior is learned via delta rule on the delivered US (not perceived pain).
%  Perceived pain = kappa * prior + (1-kappa) * sensory input.
    nTrials = size(CS, 1);
    w = [0.5; 0.5];
    V = zeros(nTrials, 1);
    for t = 1:nTrials
        x = CS(t,:)';
        mu_prior = w' * x;               % top-down prior
        s = US(t);                        % bottom-up sensory signal
        V(t) = kappa * mu_prior + (1 - kappa) * s;

        % Update weights using delivered US (not perceived pain)
        delta = US(t) - mu_prior;
        w = w + eta * x * delta;
    end
end


%% ========================================================================
function [LL, b0, b1, sigma] = rescaled_LL(V, CR)
%RESCALED_LL  Log-likelihood after linear rescaling V -> CR.
    n  = numel(CR);
    X  = [ones(n,1), V(:)];
    b  = X \ CR(:);
    b0 = b(1);  b1 = b(2);
    pred  = X * b;
    resid = CR(:) - pred;
    sigma = sqrt(mean(resid.^2));
    LL    = -0.5*n*log(2*pi) - n*log(sigma) - 0.5*sum((resid/sigma).^2);
end
