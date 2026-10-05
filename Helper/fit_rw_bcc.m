function results = fit_rw_bcc(subT)

%  Rescorla-Wagner learning model with naive initialisation (.5), fit across
%  conditioning and test phases. Inlcudes cue integration (bcc).
%      V(t)   = CS(t) * w_t
%      w_{t+1} = w_t + eta * CS(t) * (US(t) - V(t))
%      R(t)   = kappa * V(t) + (1 - kappa) * US(t)
%  Weights reset to naive (0.5, 0.5) at each block boundary.


    subT = sortrows(subT, 'TrialGlobal');
    blocks = string(subT.Block);

    CS = [subT.x_face, subT.x_house];
    US = subT.TargetVAS / 100;
    CR = subT.VASRating;
    nTrials = height(subT);

    w_0      = [0.5; 0.5]; % naive, no prior association
    resetIdx = find(blocks(2:end) ~= blocks(1:end-1)) + 1;

    % Grid search over eta
    etaGrid = linspace(0.01, 1, 50);
    kappaGrid = linspace(0, 1, 30);
    
    bestLL  = -Inf; bestEta = NaN; bestKappa = NaN;

    for k = 1:numel (kappaGrid)
        for i = 1:numel(etaGrid)
            eta = etaGrid(i);
            kappa = kappaGrid (k);
            [R, ~]   = rw_forward(CS, US, eta, kappa, w_0, resetIdx);
            LL  = rescaled_LL(R, CR);
            if LL > bestLL
                bestLL  = LL; bestEta = eta; bestKappa = kappa; 
            end
        end
    end

    % Refine with fmincon
    obj      = @(p) -rescaled_LL(rw_forward(CS, US, p(1), p(2), w_0, resetIdx), CR);
    opts_opt = optimoptions('fmincon','Display','off');
    pOpt     = fmincon(obj, [bestEta bestKappa], [],[],[],[], [0.001 0], [1 1], [], opts_opt);
    etaOpt   = pOpt(1);kappaOpt = pOpt(2);

    [R, W, V] = rw_forward(CS, US, etaOpt, kappaOpt, w_0, resetIdx);
    [LL, b0, b1, sigma] = rescaled_LL(R, CR);

    % Store results
    results.w_0    = w_0;
    results.w_t    = W;
    results.eta    = etaOpt;
    results.kappa  = kappaOpt; 
    results.R      = R;
    results.V = V;
    results.CRpred = b0 + b1 * R;
    results.LL     = LL;
    results.nPar   = 5; % eta, kappa, b0, b1, sigma
    results.BIC    = -2*LL + results.nPar * log(nTrials);
    results.b0     = b0;
    results.b1     = b1;
    results.sigma  = sigma;
end


%% ========================================================================
function [R, W, V] = rw_forward(CS, US, eta, kappa, w_0, resetIdx)
    nTrials = size(CS, 1);
    w_t = w_0;                       % initialize naive weights
    V = zeros(nTrials, 1);
    R = zeros(nTrials, 1);
    W = zeros(nTrials, 2);          % weight vector 
    for t = 1:nTrials
        if any(t == resetIdx); w_t = w_0; end      % reset at block boundary
        W(t,:) = w_t';
        x = CS(t,:)';
        s = US(t); 
        V(t) = w_t' * x;
        R(t) = kappa * V(t) + (1-kappa) * s;
        % weight updating
        delta = US(t) - V(t);
        w_t = w_t + eta * x * delta;
    end
end
