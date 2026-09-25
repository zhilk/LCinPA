nStudies = numel(behavioralData);
allStudyResults = cell(nStudies, 1);

for s = 1:nStudies
    allStudyResults{s} = run_study_models(behavioralData{s}, studyLabels{s});
end

%% ── Cross-study BIC summary ──────────────────────────────────────
fprintf('\n%s\n', repmat('=',1,90));
fprintf('MODEL COMPARISON SUMMARY (ALL STUDIES)\n');
fprintf('%s\n', repmat('=',1,90));
fprintf('%-25s %7s %7s %7s | %6s %6s %6s \n', ...
    'Study','RW','BCC','LCM', 'RW','BCC','LCM');
fprintf('%-25s %7s %7s %7s | %6s %6s %6s \n', ...
    '','───BIC sum───','','', '───% best───','','');
fprintf('%s\n', repmat('-',1,90));

for s = 1:nStudies
    R = allStudyResults{s};
    fprintf('%-25s %7.0f %7.0f %7.0f | %5.1f%% %5.1f%% %5.1f%% \n', ...
        studyLabels{s}, R.sumBIC_RW, R.sumBIC_BCC, R.sumBIC_LCM, ...
        R.pctBest_RW, R.pctBest_BCC, R.pctBest_LCM);
end
fprintf('%s\n', repmat('=',1,90));

%% ========================================================================
function R = run_study_models(T, studyLabel)

    subs = unique(T.SubID);
    nSub = numel(subs);

    fprintf('\n%s\n%s — %d subjects\n%s\n', ...
        repmat('=',1,70), studyLabel, nSub, repmat('=',1,70));

    % Preallocate per-subject results
    BIC_RW  = nan(nSub, 1);
    BIC_BCC = nan(nSub, 1);
    BIC_LCM = nan(nSub, 1);
    LL_RW   = nan(nSub, 1);
    LL_BCC  = nan(nSub, 1);
    LL_LCM  = nan(nSub, 1);
    eta_RW  = nan(nSub, 1);
    kappa_RW = nan(nSub, 1);
    kappa_BCC   = nan(nSub, 1);
    kappa_LCM  = nan(nSub, 1);
    bestModel = cell(nSub, 1);

    % Per-subject storage
    allRW  = cell(nSub, 1);
    allBCC = cell(nSub, 1);
    allLCM = cell(nSub, 1);

    fprintf('Fitting models (RW, BCC, LCM): ');
    for i = 1:nSub

        sid = subs(i);
        subT = T(T.SubID == sid, :);
        fprintf('fitting subject %d\n', sid);

        % Study 3 sub 36 excluded: <reason>
        if strcmp(studyLabel, 'Study 3 (fMRI)') && sid == 36
            fprintf('skipping subject %d (excluded)\n', sid);
            continue
        end

        % remove missing trials & non-responses 
        subT = subT(subT.VASResponse==1 & ~isnan(subT.VASRating), :);

        % Fit all models
        resRW  = fit_rw(subT);
        resBCC = fit_bcc(subT);
        resLCM = fit_lcm(subT);     % full carry-over

        % Store model comparison results
        BIC_RW(i)  = resRW.BIC;
        BIC_BCC(i) = resBCC.BIC;
        BIC_LCM(i) = resLCM.BIC;
        LL_RW(i)   = resRW.LL;
        LL_BCC(i)  = resBCC.LL;
        LL_LCM(i)  = resLCM.LL;

        % model parameter
        eta_RW(i)  = resRW.eta;
        kappa_RW(i) = resRW.kappa;
        kappa_BCC(i)   = resBCC.kappa;
        kappa_LCM(i) = resLCM.kappa;
        allRW{i}   = resRW;
        allBCC{i}  = resBCC;
        allLCM{i}  = resLCM;

        % Determine winning model
        bics = [resRW.BIC, resBCC.BIC, resLCM.BIC];
        [~, best] = min(bics);
        labels = {'RW','BCC','LCM'};
        bestModel{i} = labels{best};

        if mod(i, 10) == 0; fprintf('Done with %d \n', i); end
    end
    fprintf('done.\n');

    % ── Per-subject table ──
    fprintf('\n%-6s %9s %9s %9s | %8s %9s %10s %10s %7s | %s\n', ...
        'SubID','BIC_RW','BIC_BCC','BIC_LCM', 'eta_RW','kappa_RW','kappa_BCC','kappa_LCM','Best');
    fprintf('%s\n', repmat('-',1,99));
    for i = 1:nSub
        fprintf('%-6d %9.1f %9.1f %9.1f | %8.3f %9.3f %10.3f %10.3f %7.3f | %s\n', ...
            subs(i), BIC_RW(i), BIC_BCC(i), BIC_LCM(i),  ...
            eta_RW(i), kappa_RW(i), kappa_BCC(i), kappa_LCM(i), bestModel{i});
    end

    % ── Aggregate statistics ──
    nBest_RW  = sum(strcmp(bestModel, 'RW'));
    nBest_BCC = sum(strcmp(bestModel, 'BCC'));
    nBest_LCM = sum(strcmp(bestModel, 'LCM'));

    fprintf('\n── Model Comparison ──\n');
    fprintf('  Sum BIC:  RW = %.1f,  BCC = %.1f,  LCM = %.1f\n', ...
        sum(BIC_RW), sum(BIC_BCC), sum(BIC_LCM));
    fprintf('  Best model count:  RW = %d (%.1f%%),  BCC = %d (%.1f%%),  LCM = %d (%.1f%%)\n', ...
        nBest_RW, 100*nBest_RW/nSub, nBest_BCC, 100*nBest_BCC/nSub, ...
        nBest_LCM, 100*nBest_LCM/nSub);

    fprintf('\n── KAPPA Parameter Distributions ──\n');
    fprintf('  Mean kappa (BCC):  %.3f +/- %.3f\n', mean(kappa_BCC), std(kappa_BCC));
    fprintf('  Mean kappa (RW):  %.3f +/- %.3f\n', mean(kappa_RW), std(kappa_RW));
    fprintf('  Mean kappa (LCM):  %.3f +/- %.3f\n', mean(kappa_LCM), std(kappa_LCM));

    % ── Store in output struct ──
    R.nSub            = nSub;
    R.subs            = subs;
    R.BIC_RW          = BIC_RW;
    R.BIC_BCC         = BIC_BCC;
    R.BIC_LCM         = BIC_LCM;
    R.LL_RW           = LL_RW;
    R.LL_BCC          = LL_BCC;
    R.LL_LCM          = LL_LCM;

    R.eta_RW          = eta_RW;
    R.kappa_RW        = kappa_RW;
    R.kappa_BCC       = kappa_BCC;
    R.kappa_LCM       = kappa_LCM;
    R.bestModel       = bestModel;
    ok = ~isnan(BIC_RW);
    R.sumBIC_RW  = sum(BIC_RW(ok));
    R.sumBIC_BCC = sum(BIC_BCC(ok));
    R.sumBIC_LCM = sum(BIC_LCM(ok));
    R.pctBest_RW      = 100 * nBest_RW / nSub;
    R.pctBest_BCC     = 100 * nBest_BCC / nSub;
    R.pctBest_LCM     = 100 * nBest_LCM / nSub;
    R.allRW           = allRW;
    R.allBCC          = allBCC;
    R.allLCM          = allLCM;

    
end
