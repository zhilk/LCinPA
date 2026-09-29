nStudies = numel(behavioralData);
allStudyResults = cell(nStudies, 1);

for s = 1:nStudies
    allStudyResults{s} = run_study_models(behavioralData{s}, studyLabels{s});
end

%% ── Cross-study BIC summary ──────────────────────────────────────
fprintf('\n%s\n', repmat('=',1,104));
fprintf('MODEL COMPARISON SUMMARY (ALL STUDIES)\n');
fprintf('%s\n', repmat('=',1,104));
fprintf('%-25s %8s %8s %8s %8s | %7s %7s %7s %7s\n', ...
    'Study','NULL','RW','BCC','LCM', 'NULL','RW','BCC','LCM');
fprintf('%-25s %8s %8s %8s %8s | %7s %7s %7s %7s\n', ...
    '','───── BIC sum ─────','','','', '──── % best ────','','','');
fprintf('%s\n', repmat('-',1,104));

for s = 1:nStudies
    R = allStudyResults{s};
    fprintf('%-25s %8.0f %8.0f %8.0f %8.0f | %6.1f%% %6.1f%% %6.1f%% %6.1f%%\n', ...
        studyLabels{s}, R.sumBIC_NULL, R.sumBIC_RW, R.sumBIC_BCC, R.sumBIC_LCM, ...
        R.pctBest_NULL, R.pctBest_RW, R.pctBest_BCC, R.pctBest_LCM);
end
fprintf('%s\n', repmat('=',1,104));

%% ========================================================================
function R = run_study_models(T, studyLabel)

    subs = unique(T.SubID);
    nSub = numel(subs);

    fprintf('\n%s\n%s — %d subjects\n%s\n', ...
        repmat('=',1,70), studyLabel, nSub, repmat('=',1,70));

    % Preallocate per-subject results
    BIC_NULL  = nan(nSub, 1);
    BIC_RW    = nan(nSub, 1);
    BIC_BCC   = nan(nSub, 1);
    BIC_LCM   = nan(nSub, 1);
    LL_NULL   = nan(nSub, 1);
    LL_RW     = nan(nSub, 1);
    LL_BCC    = nan(nSub, 1);
    LL_LCM    = nan(nSub, 1);
    eta_RW    = nan(nSub, 1);
    kappa_RW  = nan(nSub, 1);
    kappa_BCC = nan(nSub, 1);
    kappa_LCM = nan(nSub, 1);
    alpha_LCM = nan(nSub, 1);
    bestModel = cell(nSub, 1);

    allNULL = cell(nSub, 1);
    allRW   = cell(nSub, 1);
    allBCC  = cell(nSub, 1);
    allLCM  = cell(nSub, 1);

    fprintf('Fitting models (NULL, RW, BCC, LCM): ');
    for i = 1:nSub

        sid  = subs(i);
        subT = T(T.SubID == sid, :);

        % Study 3 sub 36 excluded: <reason>
        if strcmp(studyLabel, 'Study 3 (fMRI)') && sid == 36
            fprintf('skipping subject %d (excluded)\n', sid);
            continue
        end

        % remove missing trials & non-responses
        subT = subT(subT.VASResponse==1 & ~isnan(subT.VASRating), :);

        % Fit all models
        resNULL = fit_null(subT);
        resRW   = fit_rw(subT);
        resBCC  = fit_bcc(subT);
        resLCM  = fit_lcm(subT);

        BIC_NULL(i) = resNULL.BIC;
        BIC_RW(i)   = resRW.BIC;
        BIC_BCC(i)  = resBCC.BIC;
        BIC_LCM(i)  = resLCM.BIC;
        LL_NULL(i)  = resNULL.LL;
        LL_RW(i)    = resRW.LL;
        LL_BCC(i)   = resBCC.LL;
        LL_LCM(i)   = resLCM.LL;

        eta_RW(i)    = resRW.eta;
        kappa_RW(i)  = resRW.kappa;
        kappa_BCC(i) = resBCC.kappa;
        kappa_LCM(i) = resLCM.kappa;
        alpha_LCM(i) = resLCM.alpha;

        allNULL{i} = resNULL;
        allRW{i}   = resRW;
        allBCC{i}  = resBCC;
        allLCM{i}  = resLCM;

        % winning model
        bics   = [resNULL.BIC, resRW.BIC, resBCC.BIC, resLCM.BIC];
        labels = {'NULL','RW','BCC','LCM'};
        [~, best] = min(bics);
        bestModel{i} = labels{best};

        if mod(i, 10) == 0; fprintf('Done with %d \n', i); end
    end
    fprintf('done.\n');

    ok   = ~isnan(BIC_RW);
    nFit = nnz(ok);

    % ── Per-subject table ──
    fprintf('\n%-6s %9s %9s %9s %9s | %8s %9s %10s %10s %10s | %s\n', ...
        'SubID','BIC_NULL','BIC_RW','BIC_BCC','BIC_LCM', ...
        'eta_RW','kappa_RW','kappa_BCC','kappa_LCM','alpha_LCM','Best');
    fprintf('%s\n', repmat('-',1,118));
    for i = 1:nSub
        if isnan(BIC_RW(i)); continue; end
        fprintf('%-6d %9.1f %9.1f %9.1f %9.1f | %8.3f %9.3f %10.3f %10.3f %10.2f | %s\n', ...
            subs(i), BIC_NULL(i), BIC_RW(i), BIC_BCC(i), BIC_LCM(i), ...
            eta_RW(i), kappa_RW(i), kappa_BCC(i), kappa_LCM(i), alpha_LCM(i), ...
            bestModel{i});
    end

    % ── Aggregate statistics ──
    nBest_NULL = sum(strcmp(bestModel(ok), 'NULL'));
    nBest_RW   = sum(strcmp(bestModel(ok), 'RW'));
    nBest_BCC  = sum(strcmp(bestModel(ok), 'BCC'));
    nBest_LCM  = sum(strcmp(bestModel(ok), 'LCM'));

    fprintf('\n── Model Comparison (n = %d) ──\n', nFit);
    fprintf('  Sum BIC:  NULL = %.1f,  RW = %.1f,  BCC = %.1f,  LCM = %.1f\n', ...
        sum(BIC_NULL(ok)), sum(BIC_RW(ok)), sum(BIC_BCC(ok)), sum(BIC_LCM(ok)));
    fprintf('  Best model count:  NULL = %d (%.1f%%),  RW = %d (%.1f%%),  BCC = %d (%.1f%%),  LCM = %d (%.1f%%)\n', ...
        nBest_NULL, 100*nBest_NULL/nFit, nBest_RW, 100*nBest_RW/nFit, ...
        nBest_BCC, 100*nBest_BCC/nFit, nBest_LCM, 100*nBest_LCM/nFit);

    % ── Does any expectation model beat the US-only baseline? ──
    dRW  = BIC_RW(ok)  - BIC_NULL(ok);
    dBCC = BIC_BCC(ok) - BIC_NULL(ok);
    dLCM = BIC_LCM(ok) - BIC_NULL(ok);
    fprintf('\n── vs NULL (negative = beats US-only) ──\n');
    fprintf('  median dBIC:    RW %+.1f,  BCC %+.1f,  LCM %+.1f\n', ...
        median(dRW), median(dBCC), median(dLCM));
    fprintf('  %% beating null: RW %.0f%%,  BCC %.0f%%,  LCM %.0f%%\n', ...
        100*mean(dRW<0), 100*mean(dBCC<0), 100*mean(dLCM<0));

    % ── store ──
    R.nSub = nSub;  R.nFit = nFit;  R.subs = subs;

    R.LL_NULL  = LL_NULL;   R.LL_RW  = LL_RW;
    R.LL_BCC   = LL_BCC;    R.LL_LCM = LL_LCM;
    R.BIC_NULL = BIC_NULL;  R.BIC_RW = BIC_RW;
    R.BIC_BCC  = BIC_BCC;   R.BIC_LCM = BIC_LCM;

    R.kappa_RW  = kappa_RW;   R.kappa_BCC = kappa_BCC;
    R.kappa_LCM = kappa_LCM;  R.eta_RW    = eta_RW;
    R.alpha_LCM = alpha_LCM;

    R.bestModel    = bestModel;
    R.sumBIC_NULL  = sum(BIC_NULL(ok));
    R.sumBIC_RW    = sum(BIC_RW(ok));
    R.sumBIC_BCC   = sum(BIC_BCC(ok));
    R.sumBIC_LCM   = sum(BIC_LCM(ok));
    R.pctBest_NULL = 100 * nBest_NULL / nFit;
    R.pctBest_RW   = 100 * nBest_RW   / nFit;
    R.pctBest_BCC  = 100 * nBest_BCC  / nFit;
    R.pctBest_LCM  = 100 * nBest_LCM  / nFit;

    R.allNULL = allNULL;  R.allRW  = allRW;
    R.allBCC  = allBCC;   R.allLCM = allLCM;
end