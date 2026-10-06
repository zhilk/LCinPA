nStudies = numel(behavioralData);
allStudyResults = cell(nStudies, 1);

for s = 1:nStudies
    allStudyResults{s} = run_study_models(behavioralData{s}, studyLabels{s});
end

%% ── Cross-study BIC summary ──────────────────────────────────────
fprintf('\n%s\n', repmat('=',1,120));
fprintf('MODEL COMPARISON SUMMARY (ALL STUDIES)\n');
fprintf('%s\n', repmat('=',1,120));
fprintf('%-25s %8s %8s %8s %8s %8s | %6s %6s %6s %6s %6s\n', ...
    'Study','NULL','RW','BCC','RW+BCC','LCM', 'NULL','RW','BCC','RW+BCC','LCM');
fprintf('%-25s %8s %8s %8s %8s %8s | %6s %6s %6s %6s %6s\n', ...
    '','─── BIC sum ───','','','','', '── % best ──','','','','');
fprintf('%s\n', repmat('-',1,120));
for s = 1:nStudies
    R = allStudyResults{s};
    fprintf('%-25s %8.0f %8.0f %8.0f %8.0f %8.0f | %5.1f%% %5.1f%% %5.1f%% %5.1f%% %5.1f%%\n', ...
        studyLabels{s}, R.sumBIC_NULL, R.sumBIC_RW, R.sumBIC_BCC, R.sumBIC_HYB, R.sumBIC_LCM, ...
        R.pctBest_NULL, R.pctBest_RW, R.pctBest_BCC, R.pctBest_HYB, R.pctBest_LCM);
end
fprintf('%s\n', repmat('=',1,120));

%% ========================================================================
function R = run_study_models(T, studyLabel)

    subs = unique(T.SubID);
    nSub = numel(subs);

    fprintf('\n%s\n%s — %d subjects\n%s\n', ...
        repmat('=',1,70), studyLabel, nSub, repmat('=',1,70));

    BIC_NULL = nan(nSub,1);  BIC_RW  = nan(nSub,1);  BIC_BCC = nan(nSub,1);
    BIC_HYB  = nan(nSub,1);  BIC_LCM = nan(nSub,1);
    LL_NULL  = nan(nSub,1);  LL_RW   = nan(nSub,1);  LL_BCC  = nan(nSub,1);
    LL_HYB   = nan(nSub,1);  LL_LCM  = nan(nSub,1);
    eta_RW   = nan(nSub,1);  kappa_BCC = nan(nSub,1);
    eta_HYB  = nan(nSub,1);  kappa_HYB = nan(nSub,1);
    alpha_LCM = nan(nSub,1);
    bestModel = cell(nSub,1);

    allNULL = cell(nSub,1); allRW = cell(nSub,1); allBCC = cell(nSub,1);
    allHYB  = cell(nSub,1); allLCM = cell(nSub,1);

    fprintf('Fitting models (NULL, RW, BCC, RW+BCC, LCM): ');
    for i = 1:nSub
        sid  = subs(i);
        subT = T(T.SubID == sid, :);

        if strcmp(studyLabel, 'Study 3 (fMRI)') && sid == 36
            fprintf('skipping subject %d (excluded)\n', sid);
            continue
        end

        subT = subT(subT.VASResponse==1 & ~isnan(subT.VASRating), :);
        subT = subT(~strcmp(subT.Phase,'conditioning'), :);

        resNULL = fit_null(subT);
        resRW   = fit_rw(subT);
        resBCC  = fit_bcc(subT);
        resHYB  = fit_rw_bcc(subT);
        resLCM  = fit_lcm(subT);

        BIC_NULL(i)=resNULL.BIC; BIC_RW(i)=resRW.BIC; BIC_BCC(i)=resBCC.BIC;
        BIC_HYB(i) =resHYB.BIC;  BIC_LCM(i)=resLCM.BIC;
        LL_NULL(i) =resNULL.LL;  LL_RW(i) =resRW.LL;  LL_BCC(i) =resBCC.LL;
        LL_HYB(i)  =resHYB.LL;   LL_LCM(i)=resLCM.LL;

        eta_RW(i)    = resRW.eta;
        kappa_BCC(i) = resBCC.kappa;
        eta_HYB(i)   = resHYB.eta;
        kappa_HYB(i) = resHYB.kappa;
        alpha_LCM(i) = resLCM.alpha;

        allNULL{i}=resNULL; allRW{i}=resRW; allBCC{i}=resBCC;
        allHYB{i} =resHYB;  allLCM{i}=resLCM;

        bics   = [resNULL.BIC, resRW.BIC, resBCC.BIC, resHYB.BIC, resLCM.BIC];
        labels = {'NULL','RW','BCC','RW+BCC','LCM'};
        [~, best] = min(bics);
        bestModel{i} = labels{best};

        if mod(i,10)==0; fprintf('Done with %d \n', i); end
    end
    fprintf('done.\n');

    ok   = ~isnan(BIC_RW);
    nFit = nnz(ok);

    % ── Per-subject table ──
    fprintf('\n%-6s %9s %9s %9s %9s %9s | %7s %9s %8s %9s %9s | %s\n', ...
        'SubID','BIC_NULL','BIC_RW','BIC_BCC','BIC_HYB','BIC_LCM', ...
        'eta_RW','kappa_BCC','eta_HYB','kappa_HYB','alpha_LCM','Best');
    fprintf('%s\n', repmat('-',1,130));
    for i = 1:nSub
        if isnan(BIC_RW(i)); continue; end
        fprintf('%-6d %9.1f %9.1f %9.1f %9.1f %9.1f | %7.3f %9.3f %8.3f %9.3f %9.2f | %s\n', ...
            subs(i), BIC_NULL(i), BIC_RW(i), BIC_BCC(i), BIC_HYB(i), BIC_LCM(i), ...
            eta_RW(i), kappa_BCC(i), eta_HYB(i), kappa_HYB(i), alpha_LCM(i), ...
            bestModel{i});
    end

    % ── Aggregate ──
    labels = {'NULL','RW','BCC','RW+BCC','LCM'};
    nB = zeros(1,5);
    for k = 1:5; nB(k) = sum(strcmp(bestModel(ok), labels{k})); end

    fprintf('\n── Model Comparison (n = %d) ──\n', nFit);
    fprintf('  Sum BIC:  NULL %.1f  RW %.1f  BCC %.1f  RW+BCC %.1f  LCM %.1f\n', ...
        sum(BIC_NULL(ok)), sum(BIC_RW(ok)), sum(BIC_BCC(ok)), ...
        sum(BIC_HYB(ok)), sum(BIC_LCM(ok)));
    fprintf('  Best count:');
    for k = 1:5; fprintf('  %s %d (%.0f%%)', labels{k}, nB(k), 100*nB(k)/nFit); end
    fprintf('\n');

    % ── vs NULL ──
    d = [BIC_RW(ok), BIC_BCC(ok), BIC_HYB(ok), BIC_LCM(ok)] - BIC_NULL(ok);
    fprintf('\n── vs NULL (negative = beats US-only) ──\n');
    fprintf('  median dBIC:    RW %+.1f  BCC %+.1f  RW+BCC %+.1f  LCM %+.1f\n', median(d,1));
    fprintf('  %% beating null: RW %.0f%%  BCC %.0f%%  RW+BCC %.0f%%  LCM %.0f%%\n', 100*mean(d<0,1));

    % ── store ──
    R.nSub = nSub;  R.nFit = nFit;  R.subs = subs;
    R.BIC_NULL=BIC_NULL; R.BIC_RW=BIC_RW; R.BIC_BCC=BIC_BCC;
    R.BIC_HYB =BIC_HYB;  R.BIC_LCM=BIC_LCM;
    R.LL_NULL =LL_NULL;  R.LL_RW =LL_RW;  R.LL_BCC =LL_BCC;
    R.LL_HYB  =LL_HYB;   R.LL_LCM=LL_LCM;
    R.eta_RW=eta_RW; R.kappa_BCC=kappa_BCC;
    R.eta_HYB=eta_HYB; R.kappa_HYB=kappa_HYB; R.alpha_LCM=alpha_LCM;
    R.bestModel = bestModel;
    R.sumBIC_NULL=sum(BIC_NULL(ok)); R.sumBIC_RW =sum(BIC_RW(ok));
    R.sumBIC_BCC =sum(BIC_BCC(ok));  R.sumBIC_HYB=sum(BIC_HYB(ok));
    R.sumBIC_LCM =sum(BIC_LCM(ok));
    R.pctBest_NULL=100*nB(1)/nFit; R.pctBest_RW =100*nB(2)/nFit;
    R.pctBest_BCC =100*nB(3)/nFit; R.pctBest_HYB=100*nB(4)/nFit;
    R.pctBest_LCM =100*nB(5)/nFit;
    R.allNULL=allNULL; R.allRW=allRW; R.allBCC=allBCC;
    R.allHYB =allHYB;  R.allLCM=allLCM;
end