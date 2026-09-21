function results = fit_bcc(subT)

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
%    subT – table for one subject, sorted by TrialGlobal, with columns:
%         x_face, x_house, TargetVAS, VASRating, Phase, Block, VisualCategory
%
%  OUTPUT:
%    results – struct with fields:
%      .kappa  – best-fit parameters
%      .V            – [nTrials x 1] model predictions (raw)
%      .CRpred       – [nTrials x 1] rescaled predictions
%      .LL, .BIC     – fit statistics
%      .nPar         – effective free parameters (eta, kappa, b0, b1, sigma)

    subT = sortrows(subT, 'TrialGlobal');

    % determine w_fixed PER BLOCK 
    blocks = subT.Block;
    if iscell(blocks); blocks = string(blocks); end
    ub = unique(blocks, 'stable');          % block labels in order

    wFixed_byblock = cell(numel(ub),1);
    for bi = 1:numel(ub)
        c = subT(blocks==ub(bi) & strcmp(subT.Phase,'conditioning'), :);
        avh = mean(c.VASRating(strcmp(c.VisualCategory,'house')), 'omitnan');
        avf = mean(c.VASRating(strcmp(c.VisualCategory,'face')),  'omitnan');
        wFixed_byblock{bi} = [avf; avh]/100;
    end

    % test trials, tagged by which block they belong to ---
    testMask = ~strcmp(subT.Phase,'conditioning'); % & subT.VASResponse==1 & subT.CatchTrial==0 & ~isnan(subT.VASRating);
    testT    = subT(testMask,:);
    tb       = testT.Block; if iscell(tb); tb = string(tb); end

    % map w_fixed to its block's 
    wFixed_seq = zeros(2, height(testT));       % 2 features x nTestTrials
    for bi = 1:numel(ub)
        wFixed_seq(:, tb==ub(bi)) = repmat(wFixed_byblock{bi}, 1, nnz(tb==ub(bi)));
    end

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
        R  = bcc_forward(CS, US, wFixed_seq, kappaGrid(j));
        LL = rescaled_LL(R, CR);
        LLgrid(j) = LL;
        if LL > bestLL
            bestLL = LL;
            bestK  = kappaGrid(j);
        end
    end

    % Refine with fmincon
    obj = @(k) -rescaled_LL(bcc_forward(CS, US, wFixed_seq, k), CR);
    opts_opt = optimoptions('fmincon','Display','off');
    kOpt = fmincon(obj, bestK, [],[],[],[], 0, 1, [], opts_opt);

    [R, V] = bcc_forward(CS, US, wFixed_seq, kOpt);
    [LL, b0, b1, sigma] = rescaled_LL(R, CR);
    
    results.wFixed = wFixed_byblock;
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

    % % ---- check: LL across the kappa grid ----
    % figure('Color','w'); plot(kappaGrid, LLgrid, '-o'); hold on
    % xline(kOpt, 'r--', sprintf('\\kappa_{opt}=%.3f', kOpt));
    % xlabel('\kappa'); ylabel('log-likelihood'); box off
    % title(sprintf('Sub %d: LL profile over \\kappa', subT.SubID(1)));
end


%% ========================================================================
function [R, V] = bcc_forward(CS, US, wFixed_seq, kappa)
    nTrials = size(CS, 1);
    V = zeros(nTrials, 1);
    R = zeros(nTrials, 1);
    for t = 1:nTrials
        s = US(t);
        V(t) = CS(t,:) * wFixed_seq(:,t);
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
