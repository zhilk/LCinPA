function results = fit_lcm(subT)

%  Latent cause model, fit by maximum likelihood across conditioning and test.
%      V(t) = sum_k post_cs(k) * mu_US(k)      (CS-only posterior, pre-update)
%      R(t) = kappa * V(t) + (1 - kappa) * US(t)
%  Causes form naively within each block and reset at block boundaries,
%  parallel to RW's naive w_0 reset.  Feature SDs fixed; alpha fit by ML.

    Kcap     = 10;
    sigma_cs = 0.2;
    sigma_us = 0.1;

    subT   = sortrows(subT, 'TrialGlobal');
    US     = subT.TargetVAS / 100;
    CS     = [subT.x_face, subT.x_house];
    CR     = subT.VASRating;
    X      = [US, CS];
    nTrials = height(subT);

    blocks   = string(subT.Block);
    ub       = unique(blocks, 'stable');
    condMask = strcmp(subT.Phase, 'conditioning');

    alphaGrid = linspace(0, 10, 50);
    kappaGrid = linspace(0, 1, 30);
    nA = numel(alphaGrid);

    % ── V for each alpha, over all trials ──
    V_all = zeros(nTrials, nA);
    for ia = 1:nA
        opts = struct('alpha',alphaGrid(ia),'stickiness',0,'a',1,'mu0',0.5, ...
                      'K',Kcap,'sigma_cs',sigma_cs,'sigma_us',sigma_us);
        V_all(:,ia) = run_blocks(X, blocks, ub, opts);
    end

    % ── joint ML over alpha and kappa ──
    bestLL = -Inf; bestIA = NaN; bestKp = NaN; LL_single = -Inf;
    for ia = 1:nA
        for kp = kappaGrid
            LL = rescaled_LL(kp*V_all(:,ia) + (1-kp)*US, CR);
            if ia==1 && LL > LL_single; LL_single = LL; end   % alpha = 0 baseline
            if LL > bestLL; bestLL = LL; bestIA = ia; bestKp = kp; end
        end
    end
    alphaOpt = alphaGrid(bestIA);
    V_best   = V_all(:,bestIA);

    % ── refine kappa at the ML alpha ──
    obj  = @(kp) -rescaled_LL(kp*V_best + (1-kp)*US, CR);
    kOpt = fmincon(obj, bestKp, [],[],[],[], 0, 1, [], ...
                   optimoptions('fmincon','Display','off'));
    R = kOpt*V_best + (1-kOpt)*US;
    [LL, b0, b1, sigma] = rescaled_LL(R, CR);

    % ── structure re-run at fixed alpha = 1 (descriptive, see methods) ──
    optsB = struct('alpha',1,'stickiness',0,'a',1,'mu0',0.5, ...
                   'K',Kcap,'sigma_cs',sigma_cs,'sigma_us',sigma_us);
    [z, Nk, SumF, condCauses] = run_blocks_struct(X, blocks, ub, condMask, optsB);
    morphMask = abs(CS(:,1)-0.5)<1e-9 & abs(CS(:,2)-0.5)<1e-9;

    % ── store ──
    results.alpha    = alphaOpt;
    results.sigma_cs = sigma_cs;      % fixed
    results.sigma_us = sigma_us;      % fixed
    results.kappa    = kOpt;
    results.V        = V_best;
    results.R        = R;
    results.CRpred   = b0 + b1*R;
    results.LL       = LL;
    results.LL_single       = LL_single;
    results.logLR_vs_single = LL - LL_single;
    results.nPar     = 5;   % alpha, kappa, b0, b1, sigma
    results.BIC      = -2*LL + results.nPar*log(nTrials);
    results.b0 = b0; results.b1 = b1; results.sigma = sigma;
    % latent structure (all trials)
    results.z          = z;
    results.Nk         = Nk;
    results.SumF       = SumF;
    results.condCauses = condCauses;
    results.morph_z    = z(morphMask);
    results.morph_newcause_frac = mean(~ismember(z(morphMask), condCauses));
    results.mu0 = 0.5; results.a = 1;
end

%% ---- V helper: one pass per block, causes fresh each block ----
function V = run_blocks(X, blocks, ub, opts)
    V = zeros(size(X,1),1);
    for bi = 1:numel(ub)
        idx = blocks==ub(bi);
        r = lcm_infer(X(idx,:), opts);
        V(idx) = r.V;
    end
end

%% ---- structure helper: per-trial z + sufficient stats, empties removed ----
function [zAll, NkAll, SumFAll, condCauses] = ...
         run_blocks_struct(X, blocks, ub, condMask, opts)
    zAll = zeros(size(X,1),1);
    NkAll = []; SumFAll = []; condCauses = []; offset = 0;
    for bi = 1:numel(ub)
        idx = find(blocks==ub(bi));
        r   = lcm_infer(X(idx,:), opts);
        act = find(r.Nk > 0);                      % used causes only
        remap = zeros(1, numel(r.Nk));
        remap(act) = 1:numel(act);
        zAll(idx) = remap(r.z)' + offset;
        ci = condMask(idx);
        for j = 1:numel(act)
            f = find(r.z == act(j), 1);            % first trial in this cause
            if ~isempty(f) && ci(f); condCauses(end+1) = j + offset; end %#ok
        end
        NkAll   = [NkAll,   r.Nk(act)];            %#ok
        SumFAll = [SumFAll; r.SumF(act,:)];        %#ok
        offset  = offset + numel(act);
    end
end
