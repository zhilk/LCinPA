function diag = diagnose_causes(allStudyResults, behavioralData, studyLabels)
%DIAGNOSE_CAUSES  Post-hoc diagnostic of LCM cause structure.
%
%   diag = diagnose_causes(allStudyResults, behavioralData, studyLabels)
%
%   Uses already-fitted LCM-reset results stored in allStudyResults
%   (from run_models) — no re-fitting required.
%
%   For each participant, extracts:
%     - Number of active causes per block
%     - What cause each mix trial gets MAP-assigned to (learning vs. novel)
%     - Whether novel-cause participants show weaker behavioral effects
%     - Posterior entropy on mix trials (assignment uncertainty)
%
%   INPUTS:
%     allStudyResults – {nStudies×1} cell of structs from run_models()
%     behavioralData  – {nStudies×1} cell of tables (raw CSV data)
%     studyLabels     – {nStudies×1} cell of study name strings
%
%   OUTPUT:
%     diag – {nStudies×1} cell of per-study diagnostic structs

nStudies = numel(allStudyResults);
diag = cell(nStudies, 1);

fprintf('\n%s\n', repmat('=', 1, 90));
fprintf('CAUSE STRUCTURE DIAGNOSTIC (block-reset LCM)\n');
fprintf('%s\n', repmat('=', 1, 90));

for s = 1:nStudies
    R = allStudyResults{s};
    T = behavioralData{s};
    subs = unique(T.SubID);
    nSub = R.nSub;

    fprintf('\n── %s (n = %d) ──\n', studyLabels{s}, nSub);

    % Per-subject storage
    nK_pla   = nan(nSub, 1);   % causes in placebo block
    nK_noc   = nan(nSub, 1);   % causes in nocebo block
    nK_total = nan(nSub, 1);   % total causes (both blocks)

    % Mix trial assignment categories
    pct_to_learning = nan(nSub, 1);  % % mix trials assigned to a learning cause
    pct_to_novel    = nan(nSub, 1);  % % mix trials assigned to a novel 3rd+ cause

    % Per-block breakdown
    pct_learning_pla = nan(nSub, 1);
    pct_learning_noc = nan(nSub, 1);

    % Posterior entropy on mix trials (higher = more uncertain assignment)
    entropy_mix = nan(nSub, 1);

    % Behavioral effect split by assignment type
    behav_effect     = R.dir_obs;
    assigned_novel   = false(nSub, 1);  % majority of mix trials to novel cause

    for i = 1:nSub
        sid = subs(i);
        subT = sortrows(T(T.SubID == sid, :), 'TrialGlobal');
        resLCM = R.allLCM_reset{i};

        if ~resLCM.resetUsed
            % Shouldn't happen for block-reset fits, but guard
            continue
        end

        nK1 = resLCM.nK_block1;
        nK2 = resLCM.nK_block2;
        nK_total(i) = nK1 + nK2;
        bc = resLCM.blockChange;

        % Determine which block is placebo vs nocebo
        blocks = subT.Block;
        if iscell(blocks); blocks = string(blocks); end

        block1_label = blocks(1);
        if block1_label == "placebo"
            nK_pla(i) = nK1;
            nK_noc(i) = nK2;
            pla_cause_ids = 1:nK1;           % cause IDs for block 1
            noc_cause_ids = (nK1+1):(nK1+nK2); % cause IDs for block 2
        else
            nK_noc(i) = nK1;
            nK_pla(i) = nK2;
            noc_cause_ids = 1:nK1;
            pla_cause_ids = (nK1+1):(nK1+nK2);
        end

        % Identify mix trials
        mixIdx = strcmp(subT.Phase, 'test') & strcmp(subT.VisualCategory, 'mix');
        if ~any(mixIdx); continue; end

        mixBlock = blocks(mixIdx);
        z_mix    = resLCM.z(mixIdx);

        % Classify mix trials within their block
        % For each block, identify the "learning causes" = causes that had
        % conditioning trials assigned to them. In a 2-cause-per-block
        % scenario, both are learning causes. A 3rd cause in a block is novel.
        %
        % Strategy: look at conditioning-phase cause assignments per block
        condIdx = strcmp(subT.Phase, 'conditioning');
        z_all   = resLCM.z;

        % Block 1 learning causes
        cond_block1 = condIdx & (1:height(subT))' < bc;
        z_cond_b1   = unique(z_all(cond_block1));

        % Block 2 learning causes
        cond_block2 = condIdx & (1:height(subT))' >= bc;
        z_cond_b2   = unique(z_all(cond_block2));

        if block1_label == "placebo"
            learning_causes_pla = z_cond_b1;
            learning_causes_noc = z_cond_b2;
        else
            learning_causes_noc = z_cond_b1;
            learning_causes_pla = z_cond_b2;
        end

        % Classify each mix trial
        nMix = numel(z_mix);
        is_learning = false(nMix, 1);

        pla_mix = mixBlock == "placebo";
        noc_mix = mixBlock == "nocebo";

        for t = 1:nMix
            if pla_mix(t)
                is_learning(t) = ismember(z_mix(t), learning_causes_pla);
            else
                is_learning(t) = ismember(z_mix(t), learning_causes_noc);
            end
        end

        pct_to_learning(i) = 100 * mean(is_learning);
        pct_to_novel(i)    = 100 * mean(~is_learning);

        if any(pla_mix)
            pct_learning_pla(i) = 100 * mean(is_learning(pla_mix));
        end
        if any(noc_mix)
            pct_learning_noc(i) = 100 * mean(is_learning(noc_mix));
        end

        % Posterior entropy on mix trials
        post_mix = resLCM.post_full(mixIdx, :);
        H = zeros(nMix, 1);
        for t = 1:nMix
            p = post_mix(t, :);
            p = p(p > 0);   % avoid log(0)
            H(t) = -sum(p .* log2(p));
        end
        entropy_mix(i) = mean(H);

        % Tag if majority assigned to novel cause
        assigned_novel(i) = pct_to_novel(i) > 50;
    end

    % ── Print cause count distribution ──
    fprintf('\n  Causes per block:\n');
    fprintf('    Placebo:  mean = %.1f, median = %.0f, range = [%.0f, %.0f]\n', ...
        mean(nK_pla, 'omitnan'), median(nK_pla, 'omitnan'), ...
        min(nK_pla), max(nK_pla));
    fprintf('    Nocebo:   mean = %.1f, median = %.0f, range = [%.0f, %.0f]\n', ...
        mean(nK_noc, 'omitnan'), median(nK_noc, 'omitnan'), ...
        min(nK_noc), max(nK_noc));
    fprintf('    Total:    mean = %.1f, median = %.0f, range = [%.0f, %.0f]\n', ...
        mean(nK_total, 'omitnan'), median(nK_total, 'omitnan'), ...
        min(nK_total), max(nK_total));

    % Histogram of per-block cause counts
    fprintf('\n  Per-block cause count distribution:\n');
    for nk = 1:5
        n_pla = sum(nK_pla == nk);
        n_noc = sum(nK_noc == nk);
        if n_pla > 0 || n_noc > 0
            fprintf('    %d causes:  placebo = %d subj (%.0f%%),  nocebo = %d subj (%.0f%%)\n', ...
                nk, n_pla, 100*n_pla/nSub, n_noc, 100*n_noc/nSub);
        end
    end

    % ── Mix trial assignment breakdown ──
    fprintf('\n  Mix trial cause assignment:\n');
    fprintf('    Assigned to a LEARNING cause:  mean = %.1f%%\n', ...
        mean(pct_to_learning, 'omitnan'));
    fprintf('    Assigned to a NOVEL cause:     mean = %.1f%%\n', ...
        mean(pct_to_novel, 'omitnan'));
    fprintf('    By block:  placebo = %.1f%% learning,  nocebo = %.1f%% learning\n', ...
        mean(pct_learning_pla, 'omitnan'), mean(pct_learning_noc, 'omitnan'));

    nNovel = sum(assigned_novel);
    fprintf('    Subjects with >50%% novel assignment:  %d/%d (%.0f%%)\n', ...
        nNovel, nSub, 100*nNovel/nSub);

    % ── Posterior entropy ──
    fprintf('\n  Mix trial posterior entropy (bits):\n');
    fprintf('    Mean = %.2f, SD = %.2f\n', ...
        mean(entropy_mix, 'omitnan'), std(entropy_mix, 'omitnan'));

    % ── Key test: does novel-cause assignment weaken the behavioral effect? ──
    fprintf('\n  Behavioral effect by assignment type:\n');
    valid = ~isnan(behav_effect);
    learn_group = valid & ~assigned_novel;
    novel_group = valid & assigned_novel;

    if sum(learn_group) >= 3 && sum(novel_group) >= 3
        mean_learn = mean(behav_effect(learn_group));
        mean_novel = mean(behav_effect(novel_group));
        fprintf('    Learning-cause group (n=%d):  mean behav effect = %+.2f\n', ...
            sum(learn_group), mean_learn);
        fprintf('    Novel-cause group   (n=%d):  mean behav effect = %+.2f\n', ...
            sum(novel_group), mean_novel);

        [p_rs, ~, stats_rs] = ranksum(behav_effect(learn_group), ...
            behav_effect(novel_group));
        if isfield(stats_rs, 'zval')
            fprintf('    Rank-sum test: z = %.2f, p = %.4f\n', stats_rs.zval, p_rs);
        else
            fprintf('    Rank-sum test: p = %.4f (exact test, n too small for z)\n', p_rs);
        end
    elseif sum(novel_group) == 0
        fprintf('    No subjects with majority novel-cause assignment.\n');
        fprintf('    All %d subjects assigned mix trials to learning causes.\n', ...
            sum(learn_group));
        p_rs = NaN;
    else
        fprintf('    Not enough subjects in both groups for comparison.\n');
        fprintf('    Learning-cause group: n=%d,  Novel-cause group: n=%d\n', ...
            sum(learn_group), sum(novel_group));
        p_rs = NaN;
    end

    % ── Correlation: entropy vs. behavioral effect ──
    fprintf('\n  Entropy vs. behavioral effect:\n');
    both_valid = ~isnan(entropy_mix) & ~isnan(behav_effect);
    rho_ent = NaN; p_ent = NaN;
    if sum(both_valid) >= 5
        [rho_ent, p_ent] = corr(entropy_mix(both_valid), ...
            behav_effect(both_valid), 'Type', 'Spearman');
        fprintf('    Spearman rho = %.3f, p = %.4f\n', rho_ent, p_ent);
        if rho_ent < 0
            fprintf('    (Negative = more uncertain assignment -> weaker effect)\n');
        else
            fprintf('    (Positive = more uncertain assignment -> stronger effect)\n');
        end
    end

    % ── Store diagnostics ──
    D.nK_pla            = nK_pla;
    D.nK_noc            = nK_noc;
    D.nK_total          = nK_total;
    D.pct_to_learning   = pct_to_learning;
    D.pct_to_novel      = pct_to_novel;
    D.pct_learning_pla  = pct_learning_pla;
    D.pct_learning_noc  = pct_learning_noc;
    D.entropy_mix       = entropy_mix;
    D.assigned_novel    = assigned_novel;
    D.behav_effect      = behav_effect;
    D.p_ranksum_groups  = p_rs;
    D.rho_entropy       = rho_ent;
    D.p_entropy         = p_ent;

    diag{s} = D;
end

fprintf('\n%s\n', repmat('=', 1, 90));
fprintf('INTERPRETATION GUIDE\n');
fprintf('%s\n', repmat('-', 1, 90));
fprintf('If most participants have 2 causes per block and mix trials are assigned\n');
fprintf('to learning causes, the cause-assignment mechanism works as described.\n\n');
fprintf('If many participants open a 3rd cause and mix trials land there,\n');
fprintf('the MAP story weakens — but the SOFT posterior can still carry\n');
fprintf('directional information because the posterior weights over the\n');
fprintf('learning causes are asymmetric (closer cause gets more weight).\n');
fprintf('%s\n', repmat('=', 1, 90));

end
