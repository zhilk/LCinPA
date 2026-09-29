function results = fit_bcc(subT)

%  Bayesian Cue Combination with fixed objective weights, fit across
%  conditioning and test phases.
%      V(t) = CS(t) * w_fixed     (w_fixed = delivered pain per CS, per block)
%      R(t) = kappa * V(t) + (1 - kappa) * US(t)

    subT   = sortrows(subT, 'TrialGlobal');
    blocks = string(subT.Block);
    cat    = string(subT.VisualCategory);
    isCond = strcmp(subT.Phase,'conditioning');
    ub     = unique(blocks,'stable');

    % fixed weights = delivered pain per CS per block
    w_byblock = zeros(2, numel(ub));
    for bi = 1:numel(ub)
        b = ub(bi);
        w_byblock(1,bi) = subT.TargetVAS(find(isCond & blocks==b & cat=="face", 1))/100;
        w_byblock(2,bi) = subT.TargetVAS(find(isCond & blocks==b & cat=="house",1))/100;
    end

    % map weights to every trial
    w_seq = zeros(2, height(subT));
    for bi = 1:numel(ub)
        m = blocks==ub(bi);
        w_seq(:, m) = repmat(w_byblock(:,bi), 1, nnz(m));
    end

    CS = [subT.x_face, subT.x_house];
    US = subT.TargetVAS / 100;
    CR = subT.VASRating;
    nTrials = height(subT);

    % grid search over kappa
    kappaGrid = linspace(0, 1, 30);
    bestLL = -Inf; bestK = NaN;
    
    for j = 1:numel(kappaGrid)
        R  = bcc_forward(CS, US, w_seq, kappaGrid(j));
        LL = rescaled_LL(R, CR);
        if LL > bestLL; bestLL = LL; bestK = kappaGrid(j); end
    end

    % refine
    obj  = @(k) -rescaled_LL(bcc_forward(CS, US, w_seq, k), CR);
    kOpt = fmincon(obj, bestK, [],[],[],[], 0, 1, [], ...
                   optimoptions('fmincon','Display','off'));

    [R, V] = bcc_forward(CS, US, w_seq, kOpt);
    [LL, b0, b1, sigma] = rescaled_LL(R, CR);

    results.w       = w_byblock;
    results.kappa   = kOpt;
    results.R       = R;
    results.V       = V;
    results.CRpred  = b0 + b1 * R;
    results.LL      = LL;
    results.nPar    = 4; % kappa, b0, b1, sigma
    results.BIC     = -2*LL + results.nPar * log(nTrials);
    results.b0 = b0; results.b1 = b1; results.sigma = sigma;
end

%% ========================================================================
function [R, V] = bcc_forward(CS, US, w_seq, kappa)
    nTrials = size(CS, 1);
    V = zeros(nTrials,1); R = zeros(nTrials,1);
    for t = 1:nTrials
        V(t) = CS(t,:) * w_seq(:,t);
        R(t) = kappa * V(t) + (1-kappa) * US(t);
    end
end

%% ========================================================================
function [LL, b0, b1, sigma] = rescaled_LL(R, CR)
    n  = numel(CR); X = [ones(n,1), R(:)];
    b  = X \ CR(:); b0 = b(1); b1 = b(2);
    resid = CR(:) - X*b;
    sigma = sqrt(mean(resid.^2));
    LL = -0.5*n*log(2*pi) - n*log(sigma) - 0.5*sum((resid/sigma).^2);
end