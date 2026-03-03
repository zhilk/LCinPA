function domCause = analyze_dominant_causes(allStudyResults, behavioralData, studyLabels)
%ANALYZE_DOMINANT_CAUSES  Extract dominant cause structure from block-reset LCM.
%
%   domCause = analyze_dominant_causes(allStudyResults, behavioralData, studyLabels)
%
%   For each participant and each block, this function:
%     1. Identifies all causes and their learned US means & CS profiles
%     2. Ranks causes by the posterior weight they receive on MIX trials
%     3. Extracts the top-2 dominant causes and their properties
%     4. Computes the posterior weight allocated to each cause type
%        (low-pain vs. high-pain) on mix trials per block
%
%   This addresses the concern that participants infer >2 causes per block
%   by showing that the two dominant causes align with the conditioning-phase
%   categories and the asymmetric assignment still holds.
%
%   INPUTS:
%     allStudyResults – {nStudies×1} cell of structs from run_models()
%     behavioralData  – {nStudies×1} cell of tables (raw CSV data)
%     studyLabels     – {nStudies×1} cell of study name strings
%
%   OUTPUT:
%     domCause – {nStudies×1} cell of per-study structs with fields:
%       .top2_us_pla      – [nSub x 2] US means of top-2 causes (placebo)
%       .top2_us_noc      – [nSub x 2] US means of top-2 causes (nocebo)
%       .top2_postw_pla   – [nSub x 2] posterior weights of top-2 (placebo)
%       .top2_postw_noc   – [nSub x 2] posterior weights of top-2 (nocebo)
%       .postw_low_pla    – [nSub x 1] posterior weight to low-pain cause (placebo mix)
%       .postw_high_pla   – [nSub x 1] posterior weight to high-pain cause (placebo mix)
%       .postw_low_noc    – [nSub x 1] posterior weight to low-pain cause (nocebo mix)
%       .postw_high_noc   – [nSub x 1] posterior weight to high-pain cause (nocebo mix)
%       .postw_other_pla  – [nSub x 1] posterior weight to non-top-2 causes (placebo)
%       .postw_other_noc  – [nSub x 1] posterior weight to non-top-2 causes (nocebo)

nStudies = numel(allStudyResults);
domCause = cell(nStudies, 1);

fprintf('\n%s\n', repmat('=', 1, 90));
fprintf('DOMINANT CAUSE ANALYSIS (block-reset LCM)\n');
fprintf('%s\n', repmat('=', 1, 90));

for s = 1:nStudies
    R = allStudyResults{s};
    T = behavioralData{s};
    subs = unique(T.SubID);
    nSub = R.nSub;

    fprintf('\n── %s (n = %d) ──\n', studyLabels{s}, nSub);

    % Preallocate
    top2_us_pla     = nan(nSub, 2);
    top2_us_noc     = nan(nSub, 2);
    top2_postw_pla  = nan(nSub, 2);
    top2_postw_noc  = nan(nSub, 2);

    postw_low_pla   = nan(nSub, 1);
    postw_high_pla  = nan(nSub, 1);
    postw_low_noc   = nan(nSub, 1);
    postw_high_noc  = nan(nSub, 1);
    postw_other_pla = nan(nSub, 1);
    postw_other_noc = nan(nSub, 1);

    for i = 1:nSub
        sid = subs(i);
        subT = sortrows(T(T.SubID == sid, :), 'TrialGlobal');
        resLCM = R.allLCM_reset{i};

        if ~resLCM.resetUsed; continue; end

        nK1 = resLCM.nK_block1;
        nK2 = resLCM.nK_block2;
        bc  = resLCM.blockChange;
        nK  = numel(resLCM.Nk);

        % Compute US mean for each cause
        Nk   = resLCM.Nk;
        SumF = resLCM.SumF;
        mu0  = resLCM.mu0;
        a_ps = resLCM.a_pseudo;

        mu_us = nan(1, nK);
        for k = 1:nK
            if Nk(k) > 0
                mu_us(k) = (SumF(k,1) + a_ps * mu0) / (Nk(k) + a_ps);
            else
                mu_us(k) = mu0;
            end
        end

        % Determine which block is first
        blocks = subT.Block;
        if iscell(blocks); blocks = string(blocks); end
        block1_label = blocks(1);

        if block1_label == "placebo"
            pla_cause_ids = 1:nK1;
            noc_cause_ids = (nK1+1):(nK1+nK2);
        else
            noc_cause_ids = 1:nK1;
            pla_cause_ids = (nK1+1):(nK1+nK2);
        end

        % Identify mix trials
        mixIdx = strcmp(subT.Phase, 'test') & strcmp(subT.VisualCategory, 'mix');
        if ~any(mixIdx); continue; end

        mixBlock = blocks(mixIdx);
        post_mix = resLCM.post_full(mixIdx, :);

        pla_mix = mixBlock == "placebo";
        noc_mix = mixBlock == "nocebo";

        % ── PLACEBO block: posterior weight on mix trials per cause ──
        if any(pla_mix) && ~isempty(pla_cause_ids)
            post_pla = post_mix(pla_mix, :);
            mean_post_pla = mean(post_pla, 1);  % [1 x nK] avg posterior per cause

            % Weight going to placebo-block causes
            w_pla_causes = mean_post_pla(pla_cause_ids);

            % Rank by posterior weight on mix trials
            [w_sorted, sort_idx] = sort(w_pla_causes, 'descend');

            % Top-2 causes
            n_top = min(2, numel(w_sorted));
            top_ids = pla_cause_ids(sort_idx(1:n_top));

            top2_us_pla(i, 1:n_top)    = mu_us(top_ids) * 100;  % back to VAS
            top2_postw_pla(i, 1:n_top) = w_sorted(1:n_top);

            % Sort top-2 by US mean (low first, high second)
            [top2_us_pla(i, 1:n_top), si] = sort(top2_us_pla(i, 1:n_top));
            top2_postw_pla(i, 1:n_top) = top2_postw_pla(i, si);

            % Posterior weight to low-pain vs high-pain cause
            % Low-pain = the one with lower US mean, high-pain = higher
            if n_top >= 2
                postw_low_pla(i)  = top2_postw_pla(i, 1);
                postw_high_pla(i) = top2_postw_pla(i, 2);
            elseif n_top == 1
                % Only one cause: assign based on US mean
                if top2_us_pla(i, 1) < 55
                    postw_low_pla(i)  = top2_postw_pla(i, 1);
                    postw_high_pla(i) = 0;
                else
                    postw_low_pla(i)  = 0;
                    postw_high_pla(i) = top2_postw_pla(i, 1);
                end
            end

            % Posterior weight to non-top-2 causes (everything else)
            postw_other_pla(i) = 1 - sum(w_sorted(1:n_top));
        end

        % ── NOCEBO block: same analysis ──
        if any(noc_mix) && ~isempty(noc_cause_ids)
            post_noc = post_mix(noc_mix, :);
            mean_post_noc = mean(post_noc, 1);

            w_noc_causes = mean_post_noc(noc_cause_ids);
            [w_sorted, sort_idx] = sort(w_noc_causes, 'descend');

            n_top = min(2, numel(w_sorted));
            top_ids = noc_cause_ids(sort_idx(1:n_top));

            top2_us_noc(i, 1:n_top)    = mu_us(top_ids) * 100;
            top2_postw_noc(i, 1:n_top) = w_sorted(1:n_top);

            [top2_us_noc(i, 1:n_top), si] = sort(top2_us_noc(i, 1:n_top));
            top2_postw_noc(i, 1:n_top) = top2_postw_noc(i, si);

            if n_top >= 2
                postw_low_noc(i)  = top2_postw_noc(i, 1);
                postw_high_noc(i) = top2_postw_noc(i, 2);
            elseif n_top == 1
                if top2_us_noc(i, 1) < 45
                    postw_low_noc(i)  = top2_postw_noc(i, 1);
                    postw_high_noc(i) = 0;
                else
                    postw_low_noc(i)  = 0;
                    postw_high_noc(i) = top2_postw_noc(i, 1);
                end
            end

            postw_other_noc(i) = 1 - sum(w_sorted(1:n_top));
        end
    end

    % ── Print summary ──
    fprintf('\n  Top-2 cause US means (VAS scale):\n');
    fprintf('    Placebo — Low cause:  M = %.1f (SD = %.1f)\n', ...
        mean(top2_us_pla(:,1), 'omitnan'), std(top2_us_pla(:,1), 'omitnan'));
    fprintf('    Placebo — High cause: M = %.1f (SD = %.1f)\n', ...
        mean(top2_us_pla(:,2), 'omitnan'), std(top2_us_pla(:,2), 'omitnan'));
    fprintf('    Nocebo  — Low cause:  M = %.1f (SD = %.1f)\n', ...
        mean(top2_us_noc(:,1), 'omitnan'), std(top2_us_noc(:,1), 'omitnan'));
    fprintf('    Nocebo  — High cause: M = %.1f (SD = %.1f)\n', ...
        mean(top2_us_noc(:,2), 'omitnan'), std(top2_us_noc(:,2), 'omitnan'));

    fprintf('\n  Posterior weight on mix trials (mean across subjects):\n');
    fprintf('    Placebo — Low-pain cause:  %.3f\n', mean(postw_low_pla, 'omitnan'));
    fprintf('    Placebo — High-pain cause: %.3f\n', mean(postw_high_pla, 'omitnan'));
    fprintf('    Placebo — Other causes:    %.3f\n', mean(postw_other_pla, 'omitnan'));
    fprintf('    Nocebo  — Low-pain cause:  %.3f\n', mean(postw_low_noc, 'omitnan'));
    fprintf('    Nocebo  — High-pain cause: %.3f\n', mean(postw_high_noc, 'omitnan'));
    fprintf('    Nocebo  — Other causes:    %.3f\n', mean(postw_other_noc, 'omitnan'));

    % Key test: in placebo, more weight to low-pain cause;
    % in nocebo, more weight to high-pain cause
    shift_pla = postw_low_pla - postw_high_pla;  % positive = more to low-pain
    shift_noc = postw_high_noc - postw_low_noc;   % positive = more to high-pain

    valid_pla = ~isnan(shift_pla);
    valid_noc = ~isnan(shift_noc);

    fprintf('\n  Asymmetry test (posterior weight difference on mix trials):\n');
    if sum(valid_pla) >= 5
        [p_pla, ~, stats_pla] = signrank(shift_pla(valid_pla));
        fprintf('    Placebo (low > high): median = %.3f, z = %.2f, p = %.4f\n', ...
            median(shift_pla(valid_pla)), stats_pla.zval, p_pla);
    end
    if sum(valid_noc) >= 5
        [p_noc, ~, stats_noc] = signrank(shift_noc(valid_noc));
        fprintf('    Nocebo  (high > low): median = %.3f, z = %.2f, p = %.4f\n', ...
            median(shift_noc(valid_noc)), stats_noc.zval, p_noc);
    end

    % ── Store ──
    D.top2_us_pla     = top2_us_pla;
    D.top2_us_noc     = top2_us_noc;
    D.top2_postw_pla  = top2_postw_pla;
    D.top2_postw_noc  = top2_postw_noc;
    D.postw_low_pla   = postw_low_pla;
    D.postw_high_pla  = postw_high_pla;
    D.postw_low_noc   = postw_low_noc;
    D.postw_high_noc  = postw_high_noc;
    D.postw_other_pla = postw_other_pla;
    D.postw_other_noc = postw_other_noc;

    domCause{s} = D;
end

fprintf('\n%s\n', repmat('=', 1, 90));
end
