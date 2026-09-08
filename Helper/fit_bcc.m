function results = fit_bcc(T)

%  Fit Bayesian Cue Combination model to one participant's data.
%
%  The BCC model (Büchel et al., 2014) computes perceived pain as a
%  precision-weighted combination of a  prior and the bottom-up sensory signal:
%
%      CR = kappa * V(t) + (1 - kappa) * s(t)
% 
%  in which the expectation V(t) is the cue vector CS with the learned
%  weights from the conditioning phase: 
%       
%      V(t) = CS(t) * w_fixed 
%
%  INPUT:
%    T – table for one subject, sorted by TrialGlobal, with columns:
%         x_face, x_house, TargetVAS, VASRating, Phase, Block, VisualCategory
%
%  OUTPUT:
%    results – struct with fields:
%      .kappa  – best-fit parameters
%      .V            – [nTrials x 1] model predictions (raw)
%      .CRpred       – [nTrials x 1] rescaled predictions
%      .LL, .BIC     – fit statistics
%      .nPar         – effective free parameters (eta, kappa, b0, b1, sigma)

    T = sortrows(T, 'TrialGlobal');

    % determine w_fixed 
    keep = strcmp(T.Phase,'conditioning') & T.VASResponse==1 & T.CatchTrial==0;
    C = T(keep,:);
    av_house = mean(C.VASRating(strcmp(C.VisualCategory,'house')), 'omitnan');
    av_face  = mean(C.VASRating(strcmp(C.VisualCategory,'face')),  'omitnan');
    w_fixed = [av_face; av_house]/100;

    % extract experimental data from test phase 
    testMask = ~strcmp(T.Phase,'conditioning'); % & T.VASResponse==1 & T.CatchTrial==0 & ~isnan(T.VASRating);
    testT = T(testMask,:);

    CS = [testT.x_face, testT.x_house];
    US = testT.TargetVAS / 100;
    CR = testT.VASRating;
    nTrials = height(testT);

    % Grid search over eta and kappa
    kappaGrid = linspace(0, 1, 30);
    bestLL    = -Inf;     
    LLgrid = zeros(size(kappaGrid));
    bestK     = NaN; % best kappa

    for j = 1:numel(kappaGrid)
        R  = bcc_forward(CS, US, w_fixed,  kappaGrid(j));
        LL = rescaled_LL(R, CR);
        LLgrid(j) = LL;
        if LL > bestLL
            bestLL = LL;
            bestK  = kappaGrid(j);
        end
    end

    % Refine with fmincon
    obj = @(k) -rescaled_LL(bcc_forward(CS, US, w_fixed, k), CR);
    opts_opt = optimoptions('fmincon','Display','off');
    kOpt = fmincon(obj, bestK, [],[],[],[], 0, 1, [], opts_opt);

    [R, V] = bcc_forward(CS, US, w_fixed, kOpt);
    [LL, b0, b1, sigma] = rescaled_LL(R, CR);
    

    results.kappa  = kOpt;
    results.R      = R; 
    results.V      = V;
    results.CRpred = b0 + b1 * R;
    results.LL     = LL;
    results.nPar   = 4;                    % kappa, b0, b1, sigma
    results.BIC    = -2*LL + results.nPar * log(nTrials);
    results.b0     = b0;
    results.b1     = b1;
    results.sigma  = sigma;

    % ---- profile check: LL across the kappa grid ----
    figure('Color','w'); plot(kappaGrid, LLgrid, '-o'); hold on
    xline(kOpt, 'r--', sprintf('\\kappa_{opt}=%.3f', kOpt));
    xlabel('\kappa'); ylabel('log-likelihood'); box off
    title(sprintf('Sub %d: LL profile over \\kappa', T.SubID(1)));
end


%% ========================================================================
function [R, V] = bcc_forward(CS, US, w_fixed, kappa)
    nTrials = size(CS, 1);
    V = zeros(nTrials, 1);
    R = zeros(nTrials, 1);
    for t = 1:nTrials
        x = CS(t,:);
        s = US(t);                        % bottom-up sensory signal
        V(t) = x * w_fixed;
        R(t) = kappa * V(t) + (1-kappa) * s;
    end
end


%% ========================================================================
function [LL, b0, b1, sigma] = rescaled_LL(R, CR)
    n  = numel(CR);
    X  = [ones(n,1), R(:)];
    b  = X \ CR(:);
    b0 = b(1);  b1 = b(2);
    pred  = X * b;
    resid = CR(:) - pred;
    sigma = sqrt(mean(resid.^2));
    LL    = -0.5*n*log(2*pi) - n*log(sigma) - 0.5*sum((resid/sigma).^2);
end
