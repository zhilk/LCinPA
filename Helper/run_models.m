function allStudyResults = run_models(behavioralData, studyLabels)
%RUN_MODELS  Fit RW, BCC, and LCM to all participants across studies.
%
%  Fits three computational models per subject:
%    1. Rescorla-Wagner (RW) – associative learning baseline
%    2. Bayesian Cue Combination (BCC) – predictive coding account
%    3. Latent Cause Model (LCM) – latent cause inference
%
%  Reports: BIC comparison, logBF (multiple causes), cause assignment
%  analysis, individual-differences correlations, and sensitivity analyses.
%
%  INPUTS:
%    behavioralData – 3x1 cell array of tables (one per study)
%    studyLabels    – 3x1 cell array of study name strings
%
%  OUTPUT:
%    allStudyResults – 3x1 cell array of per-study result structs

nStudies = numel(behavioralData);
allStudyResults = cell(nStudies, 1);

for s = 1:nStudies
    allStudyResults{s} = run_study_models(behavioralData{s}, studyLabels{s});
end

%% ── Cross-study BIC summary ──────────────────────────────────────
fprintf('\n%s\n', repmat('=',1,90));
fprintf('MODEL COMPARISON SUMMARY (ALL STUDIES)\n');
fprintf('%s\n', repmat('=',1,90));
fprintf('%-25s %7s %7s %7s | %6s %6s %6s | %7s\n', ...
    'Study','RW','BCC','LCM', 'RW','BCC','LCM', 'logBF');
fprintf('%-25s %7s %7s %7s | %6s %6s %6s | %7s\n', ...
    '','───BIC sum───','','', '───% best───','','','(LCMvRW)');
fprintf('%s\n', repmat('-',1,90));

for s = 1:nStudies
    R = allStudyResults{s};
    fprintf('%-25s %7.0f %7.0f %7.0f | %5.1f%% %5.1f%% %5.1f%% | %7.1f\n', ...
        studyLabels{s}, R.sumBIC_RW, R.sumBIC_BCC, R.sumBIC_LCM, ...
        R.pctBest_RW, R.pctBest_BCC, R.pctBest_LCM, ...
        mean(R.logBF_LCM_vs_RW));
end
fprintf('%s\n', repmat('=',1,90));

%% ── Cross-study cause assignment summary ─────────────────────────
fprintf('\nCAUSE ASSIGNMENT ANALYSIS — block-reset LCM (ALL STUDIES)\n');
fprintf('%s\n', repmat('-',1,85));
fprintf('%-25s | %10s %6s | %8s %8s\n', ...
    'Study', 'MeanDiff','%Corr', 'Wil. p','Corr rho');
fprintf('%s\n', repmat('-',1,85));
for s = 1:nStudies
    R = allStudyResults{s};
    fprintf('%-25s | %+10.6f %5.1f%% | %8.4f %+8.3f\n', ...
        studyLabels{s}, ...
        R.dir_assign_reset_mean, R.dir_assign_reset_pct, ...
        R.pval_wsr_assign, R.rho_soft);
end
fprintf('%s\n', repmat('-',1,85));
fprintf('MeanDiff = MAP cause US mean (nocebo - placebo); positive = correct.\n');
fprintf('Corr rho = Spearman correlation between soft posterior shift and observed effect.\n');
fprintf('%s\n', repmat('=',1,85));

end


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
    logBF   = nan(nSub, 1);
    eta_RW  = nan(nSub, 1);
    kappa_RW = nan(nSub, 1);
    kappa_BCC   = nan(nSub, 1);
    alpha   = nan(nSub, 1);
    bestModel = cell(nSub, 1);

    % Cause assignment analysis (block-reset LCM)
    dir_obs          = nan(nSub, 1);  % observed VAS diff (noc - pla)
    dir_assign_reset = nan(nSub, 1);  % hard MAP cause US mean diff
    dir_soft_reset   = nan(nSub, 1);  % soft posterior-weighted US mean diff

    % Per-subject storage
    allRW  = cell(nSub, 1);
    allBCC = cell(nSub, 1);
    allLCM = cell(nSub, 1);
    allLCM_reset = cell(nSub, 1);

    fprintf('Fitting models (RW, BCC, LCM, LCM-reset): ');
    for i = 1:nSub
        sid = subs(i);
        subT = T(T.SubID == sid, :);

        % remove missing trials & non-responses 
        subT = subT(subT.VASResponse==1 & ~isnan(subT.VASRating), :);

        % Fit all models
        resRW  = fit_rw(subT);
        resBCC = fit_bcc(subT);
        resLCM = fit_lcm(subT);     % full carry-over
        resLCM_reset = fit_lcm(subT, true); % block reset

        % Store model comparison results
        BIC_RW(i)  = resRW.BIC;
        BIC_BCC(i) = resBCC.BIC;
        BIC_LCM(i) = resLCM.BIC;
        LL_RW(i)   = resRW.LL;
        LL_BCC(i)  = resBCC.LL;
        LL_LCM(i)  = resLCM.LL;
        logBF(i)   = resLCM.logBF_vs_rw;

        % model parameter
        eta_RW(i)  = resRW.eta;
        kappa_RW(i) = resRW.kappa;
        kappa_BCC(i)   = resBCC.kappa;
        alpha(i)   = resLCM.alpha_map;
        allRW{i}   = resRW;
        allBCC{i}  = resBCC;
        allLCM{i}  = resLCM;
        allLCM_reset{i} = resLCM_reset;

        % Determine winning model
        bics = [resRW.BIC, resBCC.BIC, resLCM.BIC];
        [~, best] = min(bics);
        labels = {'RW','BCC','LCM'};
        bestModel{i} = labels{best};

        % ── Cause assignment analysis (block-reset LCM) ──
        subT_sorted = sortrows(subT, 'TrialGlobal');
        mixIdx = strcmp(subT_sorted.Phase, 'test') & ...
                 strcmp(subT_sorted.VisualCategory, 'mix');

        if any(mixIdx)
            mixBlock = subT_sorted.Block(mixIdx);
            if iscell(mixBlock); mixBlock = string(mixBlock); end

            mixVAS = subT_sorted.VASRating(mixIdx);

            noc_idx = mixBlock == "nocebo";
            pla_idx = mixBlock == "placebo";

            % Observed behavioral effect
            if any(noc_idx) && any(pla_idx)
                dir_obs(i) = mean(mixVAS(noc_idx)) - mean(mixVAS(pla_idx));
            end

            % Compute US mean for each cause (block-reset)
            z_mix = resLCM_reset.z(mixIdx);
            Nk    = resLCM_reset.Nk;
            SumF  = resLCM_reset.SumF;
            mu0   = resLCM_reset.mu0;
            a_ps  = resLCM_reset.a_pseudo;

            nCk = numel(Nk);
            mu_us = nan(1, nCk);
            for k = 1:nCk
                if Nk(k) > 0
                    mu_us(k) = (SumF(k,1) + a_ps * mu0) / (Nk(k) + a_ps);
                else
                    mu_us(k) = mu0;
                end
            end

            % Hard MAP: look up assigned cause US mean per trial
            assigned_us = nan(numel(z_mix), 1);
            for t = 1:numel(z_mix)
                zk = z_mix(t);
                if zk <= nCk
                    assigned_us(t) = mu_us(zk);
                end
            end

            if any(noc_idx) && any(pla_idx) && ...
               any(~isnan(assigned_us(noc_idx))) && any(~isnan(assigned_us(pla_idx)))
                dir_assign_reset(i) = mean(assigned_us(noc_idx), 'omitnan') ...
                                    - mean(assigned_us(pla_idx), 'omitnan');
            end

            % Soft posterior-weighted US mean per trial
            post_mix = resLCM_reset.post_full(mixIdx, :);
            soft_us = nan(numel(z_mix), 1);
            for t = 1:numel(z_mix)
                soft_us(t) = post_mix(t, :) * mu_us(:);
            end

            if any(noc_idx) && any(pla_idx)
                dir_soft_reset(i) = mean(soft_us(noc_idx)) - mean(soft_us(pla_idx));
            end
        end

        if mod(i, 10) == 0; fprintf('%d ', i); end
    end
    fprintf('done.\n');

        % ── Per-subject table ──
    fprintf('\n%-6s %8s %8s %8s | %7s %7s %8s %8s %6s | %s\n', ...
        'SubID','BIC_RW','BIC_BCC','BIC_LCM','logBF','eta RW','kappa RW','kappa BCC','alpha','Best');
    fprintf('%s\n', repmat('-',1,84));
    for i = 1:nSub
        fprintf('%-6d %8.1f %8.1f %8.1f | %7.2f %7.3f %8.3f %8.3f %6.3f | %s\n', ...
            subs(i), BIC_RW(i), BIC_BCC(i), BIC_LCM(i), logBF(i), ...
            eta_RW(i), kappa_RW(i), kappa_BCC(i), alpha(i), bestModel{i});
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

    fprintf('\n── Evidence for Multiple Latent Causes ──\n');
    fprintf('  Mean logBF (LCM vs single-cause):  %.2f\n', mean(logBF));
    fprintf('  %% subjects with logBF > 3 (strong evidence):  %.1f%%\n', 100*mean(logBF > 3));
    fprintf('  Mean alpha:  %.3f +/- %.3f\n', mean(alpha), std(alpha));

    fprintf('\n── Parameter Distributions ──\n');
    fprintf('  Mean kappa (BCC):  %.3f +/- %.3f\n', mean(kappa_BCC), std(kappa_BCC));
    fprintf('  Mean eta (RW):  %.3f +/- %.3f\n', mean(eta_RW), std(eta_RW));

    % ── Cause assignment analysis (block-reset LCM) ──
    fprintf('\n── Cause Assignment Analysis (block-reset LCM) ──\n');
    fprintf('  With block reset, each block has independent causes.\n');
    fprintf('  VAS 50 is assigned to the low-pain cause in placebo (closer to 40)\n');
    fprintf('  and the high-pain cause in nocebo (closer to 60).\n\n');

    % Hard MAP
    fprintf('  MAP cause assignment (nocebo - placebo US mean):\n');
    n_correct_assign = sum(dir_assign_reset > 0 & ~isnan(dir_assign_reset));
    n_valid_assign = sum(~isnan(dir_assign_reset));
    fprintf('    Mean diff:           %+.6f\n', mean(dir_assign_reset, 'omitnan'));
    fprintf('    %% subjects correct: %.1f%%  (%d/%d)\n', ...
        100*n_correct_assign/max(n_valid_assign,1), n_correct_assign, n_valid_assign);

    % Wilcoxon signed-rank test
    valid_assign = dir_assign_reset(~isnan(dir_assign_reset));
    pval_wsr = NaN;
    zval_wsr = NaN;
    if numel(valid_assign) >= 5
        [pval_wsr, ~, stats_wsr] = signrank(valid_assign);
        zval_wsr = stats_wsr.zval;
        fprintf('    Wilcoxon signed-rank: z = %.2f, p = %.4f\n', zval_wsr, pval_wsr);
    end

    % Individual differences: soft posterior correlates with behavior
    fprintf('\n  Individual differences (soft posterior vs. observed effect):\n');
    rho_soft = NaN;
    pval_corr_soft = NaN;
    both_valid_soft = ~isnan(dir_soft_reset) & ~isnan(dir_obs);
    if sum(both_valid_soft) >= 5
        [rho_soft, pval_corr_soft] = corr(dir_soft_reset(both_valid_soft), ...
            dir_obs(both_valid_soft), 'Type', 'Spearman');
        fprintf('    Spearman rho = %.3f, p = %.4f\n', rho_soft, pval_corr_soft);
    end

    % Control correlations: do simpler model parameters predict the behavioral effect?
    fprintf('\n  Control correlations (simpler predictors vs. observed effect):\n');
    rho_eta = NaN; pval_eta = NaN;
    rho_kappa = NaN; pval_kappa = NaN;
    rho_alpha = NaN; pval_alpha = NaN;

    valid_obs = ~isnan(dir_obs);

    % RW learning rate vs. behavioral effect
    both_valid = valid_obs & ~isnan(eta_RW);
    if sum(both_valid) >= 5
        [rho_eta, pval_eta] = corr(eta_RW(both_valid), ...
            dir_obs(both_valid), 'Type', 'Spearman');
        fprintf('    RW eta   vs. behav: rho = %.3f, p = %.4f\n', rho_eta, pval_eta);
    end

    % BCC kappa vs. behavioral effect
    both_valid = valid_obs & ~isnan(kappa_BCC);
    if sum(both_valid) >= 5
        [rho_kappa, pval_kappa] = corr(kappa_BCC(both_valid), ...
            dir_obs(both_valid), 'Type', 'Spearman');
        fprintf('    BCC kappa vs. behav: rho = %.3f, p = %.4f\n', rho_kappa, pval_kappa);
    end

    % LCM alpha vs. behavioral effect
    both_valid = valid_obs & ~isnan(alpha);
    if sum(both_valid) >= 5
        [rho_alpha, pval_alpha] = corr(alpha(both_valid), ...
            dir_obs(both_valid), 'Type', 'Spearman');
        fprintf('    LCM alpha vs. behav: rho = %.3f, p = %.4f\n', rho_alpha, pval_alpha);
    end

    % ── Sensitivity: block reset BIC ──
    fprintf('\n── Sensitivity: LCM block reset ──\n');
    BIC_LCM_reset = nan(nSub, 1);
    for i = 1:nSub
        BIC_LCM_reset(i) = allLCM_reset{i}.BIC;
    end
    nBetter_carry = sum(BIC_LCM < BIC_LCM_reset);
    nBetter_reset = sum(BIC_LCM_reset < BIC_LCM);
    fprintf('  Carry-over better: %d/%d subjects (%.1f%%)\n', ...
        nBetter_carry, nSub, 100*nBetter_carry/nSub);
    fprintf('  Block-reset better: %d/%d subjects (%.1f%%)\n', ...
        nBetter_reset, nSub, 100*nBetter_reset/nSub);
    fprintf('  Mean BIC diff (reset - carry-over): %.2f\n', ...
        mean(BIC_LCM_reset - BIC_LCM));

    % ── Sensitivity: temporal decay ──
    fprintf('\n── Sensitivity: LCM with temporal decay ──\n');
    BIC_LCM_decay  = nan(nSub, 1);
    lambda_decay    = nan(nSub, 1);
    fprintf('  Fitting: ');
    for i = 1:nSub
        sid = subs(i);
        subT = T(T.SubID == sid, :);
        resDecay = fit_lcm_decay(subT, false);
        BIC_LCM_decay(i) = resDecay.BIC;
        lambda_decay(i)   = resDecay.lambda;
        if mod(i, 10) == 0; fprintf('%d ', i); end
    end
    fprintf('done.\n');
    nBetter_nodecay = sum(BIC_LCM < BIC_LCM_decay);
    nBetter_decay   = sum(BIC_LCM_decay < BIC_LCM);
    fprintf('  Standard LCM better: %d/%d subjects (%.1f%%)\n', ...
        nBetter_nodecay, nSub, 100*nBetter_nodecay/nSub);
    fprintf('  Decay LCM better:    %d/%d subjects (%.1f%%)\n', ...
        nBetter_decay, nSub, 100*nBetter_decay/nSub);
    fprintf('  Mean BIC diff (decay - standard): %.2f\n', ...
        mean(BIC_LCM_decay - BIC_LCM));
    fprintf('  Mean best lambda: %.3f +/- %.3f\n', mean(lambda_decay), std(lambda_decay));

    % ── Store in output struct ──
    R.nSub            = nSub;
    R.subs            = subs;
    R.BIC_RW          = BIC_RW;
    R.BIC_BCC         = BIC_BCC;
    R.BIC_LCM         = BIC_LCM;
    R.BIC_LCM_reset   = BIC_LCM_reset;
    R.BIC_LCM_decay   = BIC_LCM_decay;
    R.lambda_decay    = lambda_decay;
    R.LL_RW           = LL_RW;
    R.LL_BCC          = LL_BCC;
    R.LL_LCM          = LL_LCM;
    R.logBF_LCM_vs_RW = logBF;

    R.eta_RW          = eta_RW;
    R.kappa_RW        = kappa_RW;
    R.kappa_BCC       = kappa_BCC;
    R.alpha           = alpha;
    R.bestModel       = bestModel;
    R.sumBIC_RW       = sum(BIC_RW);
    R.sumBIC_BCC      = sum(BIC_BCC);
    R.sumBIC_LCM      = sum(BIC_LCM);
    R.pctBest_RW      = 100 * nBest_RW / nSub;
    R.pctBest_BCC     = 100 * nBest_BCC / nSub;
    R.pctBest_LCM     = 100 * nBest_LCM / nSub;
    R.allRW           = allRW;
    R.allBCC          = allBCC;
    R.allLCM          = allLCM;
    R.allLCM_reset    = allLCM_reset;

    % Cause assignment results
    R.dir_obs              = dir_obs;
    R.dir_assign_reset     = dir_assign_reset;
    R.dir_soft_reset       = dir_soft_reset;
    R.dir_assign_reset_mean = mean(dir_assign_reset, 'omitnan');
    R.dir_assign_reset_pct  = 100 * mean(dir_assign_reset > 0, 'omitnan');
    R.pval_wsr_assign      = pval_wsr;
    R.zval_wsr_assign      = zval_wsr;
    R.rho_soft             = rho_soft;
    R.pval_corr_soft       = pval_corr_soft;

    % Control correlations
    R.rho_eta              = rho_eta;
    R.pval_eta             = pval_eta;
    R.rho_kappa            = rho_kappa;
    R.pval_kappa           = pval_kappa;
    R.rho_alpha            = rho_alpha;
    R.pval_alpha           = pval_alpha;
end
