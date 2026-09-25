function results = fit_lcm(subT)
%  LCM fit: alpha, sigma_cs, sigma_us, kappa by ML.
%  CR = kappa*V + (1-kappa)*US.
%  Conditioning forms causes; at the test boundary cause COUNTS are hard-reset
%  to n_prior (means preserved) so conditioning acts as a prior — parallel to
%  how RW/BCC use conditioning only to set w0.  LL scored on TEST trials.

    Kcap    = 10;
    n_prior = 1;     % prior strength of conditioning causes entering test
                     % (fixed; matches RW/BCC single-prior treatment)

    subT = sortrows(subT, 'TrialGlobal');
    US = subT.TargetVAS / 100;
    CS = [subT.x_face, subT.x_house];
    CR = subT.VASRating;
    X  = [US, CS];

    blocks = subT.Block; if iscell(blocks); blocks = string(blocks); end
    ub = unique(blocks, 'stable');
    condMask = strcmp(subT.Phase, 'conditioning');
    testMask = ~condMask;

    US_test = US(testMask);
    CR_test = CR(testMask);
    nTest   = nnz(testMask);

    sigma_cs = 0.2;
    sigma_us = 0.1;
    alphaGrid = linspace(0, 5, 50);
    kappaGrid = linspace(0, 1, 30);
    nA = numel(alphaGrid);

    bestMargLL = -Inf; best = struct(); LL_single = -Inf;

    % V for each alpha (single sigma now)
    Vt_all = zeros(nTest, nA);
    for ia = 1:nA
        opts = struct('alpha',alphaGrid(ia),'stickiness',0,'a',1,'mu0',0.5,...
                      'K',Kcap,'sigma_cs',sigma_cs,'sigma_us',sigma_us);
        Vt = run_blocks(X, blocks, ub, condMask, opts, n_prior);
        Vt_all(:,ia) = Vt(testMask);
    end
    for kp = kappaGrid
        LLa = zeros(1,nA);
        for ia = 1:nA
            LLa(ia) = rescaled_LL(kp*Vt_all(:,ia) + (1-kp)*US_test, CR_test);
        end
        if LLa(1) > LL_single; LL_single = LLa(1); end
        margLL = logsumexp(LLa) - log(nA);
        if margLL > bestMargLL
            bestMargLL = margLL;
            best = struct('kp',kp,'Vt_all',{Vt_all},'LLa',LLa);
        end
    end

    % ── alpha-marginalized V at best (sigma, kappa) ──
    logPost = best.LLa - logsumexp(best.LLa);
    postA   = exp(logPost(:));                 % posterior over alpha
    Vt_avg  = best.Vt_all * postA;             % [nTest x 1] marginalized V

    % ── refine kappa on marginalized V ──
    obj  = @(kp) -rescaled_LL(kp*Vt_avg + (1-kp)*US_test, CR_test);
    kOpt = fmincon(obj, best.kp, [],[],[],[], 0, 1, [], ...
                   optimoptions('fmincon','Display','off'));
    R = kOpt*Vt_avg + (1-kOpt)*US_test;
    [LL, b0, b1, sigma] = rescaled_LL(R, CR_test);

    % ── modal alpha for the structure re-run (point estimate, structure only) ──
    [~, iA]     = max(best.LLa);
    alpha_mode  = alphaGrid(iA);

    % ── structure re-run: modal alpha, FIXED sigma_cs=0.25 (fitted sigma shatters) ──
    optsB = struct('alpha',1,'stickiness',0,'a',1,'mu0',0.5, ...
                   'K',Kcap,'sigma_cs',sigma_cs,'sigma_us',sigma_us);

    % ... (alpha=1 fixed for structure, as we established, since modal alpha->0 collapses causes)
    [zTest, postTest, condCauses, NkAll, SumFAll] = ...
        run_blocks_struct(X, blocks, ub, condMask, optsB, n_prior);
    CS_test = CS(testMask,:);
    morphT  = abs(CS_test(:,1)-0.5)<1e-9 & abs(CS_test(:,2)-0.5)<1e-9;

    % clear empty causes
    act = NkAll > 0;
    NkAll    = NkAll(act);
    SumFAll  = SumFAll(act,:);
    postTest = postTest(:, act);              % remap post columns too

    oldIds = find(act);
    remap  = zeros(1, numel(act)); remap(oldIds) = 1:numel(oldIds);
    zTest      = remap(zTest); zTest = zTest(:);   % column
    condCauses = remap(condCauses);
    condCauses = condCauses(condCauses>0);

    % ── store ──
    results.alpha_mode = alpha_mode;   % modal alpha (structure only; NOT a fit)
    results.alphaPost  = postA(:)';    % posterior over alphaGrid
    results.alphaGrid  = alphaGrid;
    results.sigma_cs = sigma_cs;   
    results.sigma_us = sigma_us;   
    results.nPar = 4;   % kappa, b0, b1, sigma  (alpha marginalized; sigmas fixed)
    results.BIC  = -2*bestMargLL + results.nPar*log(nTest);
    results.kappa      = kOpt;
    results.n_prior    = n_prior;
    results.Nk         = NkAll;
    results.SumF       = SumFAll;
    results.V          = Vt_avg;        % marginalized V
    results.R          = R;
    results.CRpred     = b0 + b1*R;
    results.LL         = bestMargLL;                  % MARGINAL LL (primary)
    results.LL_point   = LL;                          % LL at refined kappa (diagnostic)
    results.LL_single  = LL_single;                   % alpha=0 baseline
    results.logBF_vs_single = bestMargLL - LL_single; % marginal Bayes-factor-ish
    results.nPar       = 4;   % sigma_cs, sigma_us, kappa, b0, b1, sigma  (alpha marginalized)
    results.BIC        = -2*bestMargLL + results.nPar*log(nTest);
    results.b0 = b0; results.b1 = b1; results.sigma = sigma;
    results.z_test     = zTest;
    results.post_test  = postTest;
    results.morph_z    = zTest(morphT);
    results.morph_newcause_frac = mean(~ismember(zTest(morphT), condCauses));
    results.condCauses = condCauses;
    results.mu0 = 0.5; results.a = 1;
end

%% ---- V helper: cond pass -> hard reset -> test pass, per block ----
function V = run_blocks(X, blocks, ub, condMask, opts, n_prior)
    V = zeros(size(X,1),1);
    for bi = 1:numel(ub)
        idx  = find(blocks==ub(bi));
        cidx = idx(condMask(idx));
        tidx = idx(~condMask(idx));

        rc = lcm_infer(X(cidx,:), opts);          % form causes
        act = rc.Nk > 0;                          % active (conditioning) causes
        [Nk0, SumF0] = reset_prior(rc.Nk, rc.SumF, n_prior);

        optsT = opts;
        optsT.initNk   = Nk0(act);                % seed ONLY active causes
        optsT.initSumF = SumF0(act,:);
        rt = lcm_infer(X(tidx,:), optsT);          % test pass, seeded
        V(tidx) = rt.V;
    end
end

%% ---- structure helper: same split, returns TEST-trial z/post + suff stats ----
function [zTest, postTest, condCauses, NkAll, SumFAll] = ...
         run_blocks_struct(X, blocks, ub, condMask, opts, n_prior)
    zTest = []; postCell = {}; condCauses = [];
    NkAll = []; SumFAll = []; offset = 0;
    for bi = 1:numel(ub)
        idx  = find(blocks==ub(bi));
        cidx = idx(condMask(idx));
        tidx = idx(~condMask(idx));

        rc = lcm_infer(X(cidx,:), opts);
        act = rc.Nk > 0; nSeed = nnz(act);
        [Nk0, SumF0] = reset_prior(rc.Nk, rc.SumF, n_prior);

        optsT = opts;
        optsT.initNk   = Nk0(act);
        optsT.initSumF = SumF0(act,:);
        rt = lcm_infer(X(tidx,:), optsT);

        Kb = size(rt.post,2);
        zTest    = [zTest; rt.z + offset];         %#ok
        postCell{end+1} = rt.post;                 %#ok
        for k = 1:nSeed; condCauses(end+1) = k + offset; end %#ok
        NkAll   = [NkAll,   rt.Nk(1:Kb)];          %#ok
        SumFAll = [SumFAll; rt.SumF(1:Kb,:)];      %#ok
        offset  = offset + Kb;
    end
    postTest = blkdiag(postCell{:});
end

%% ---- hard reset: preserve mean, drop count to n_prior ----
function [Nk, SumF] = reset_prior(Nk, SumF, n_prior)
    for k = find(Nk > 0)
        s = n_prior / Nk(k);
        SumF(k,:) = SumF(k,:) * s;   % keeps SumF/Nk unchanged
        Nk(k)     = n_prior;
    end
end

%% ---- OLS rescaling likelihood ----
function [LL, b0, b1, sigma] = rescaled_LL(V, CR)
    n = numel(CR); Xr = [ones(n,1), V(:)];
    b = Xr \ CR(:); b0 = b(1); b1 = b(2);
    resid = CR(:) - Xr*b;
    sigma = sqrt(mean(resid.^2));
    LL = -0.5*n*log(2*pi) - n*log(sigma) - 0.5*sum((resid/sigma).^2);
end

%% logsum function 
function s = logsumexp(x)
    m = max(x); s = m + log(sum(exp(x - m)));
end