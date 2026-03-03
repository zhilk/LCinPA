function plot_modeling(allStudyResults, studyLabels, causeDiag)
%PLOT_MODELING  Publication-quality figures for LCPA computational modeling.
%
%   plot_modeling(allStudyResults, studyLabels)
%   plot_modeling(allStudyResults, studyLabels, causeDiag)
%
%   INPUTS:
%     allStudyResults – {3×1} cell array of structs from run_models()
%     studyLabels     – {3×1} cell array of study name strings
%     causeDiag       – (optional) {3×1} cell array from diagnose_causes()
%                       If provided, enables Panels E & F with proper
%                       learning-vs-novel group splits and entropy data.
%
%   Produces two figures:
%     Figure 1 – Model Comparison (BIC, logBF, block-reset, temporal decay)
%     Figure 2 – Cause Assignment Mechanism (5 panels, PNAS format)
%                A: schematic, B: individual diffs, C: control correlations,
%                D: learning-vs-novel group split, E: entropy correlation
%
%   Requires: al_goodplot.m on the MATLAB path.

    if nargin < 3; causeDiag = []; end

    nStudies = numel(allStudyResults);

    %% ---- Colour definitions ----
    % Model identity colours (distinct from behavioural blue/red)
    col_RW   = [0.85, 0.65, 0.13];   % amber
    col_BCC  = [0.17, 0.63, 0.53];   % teal
    col_LCM  = [0.55, 0.30, 0.65];   % purple

    % Study-graded purples (for cause-assignment panels)
    study_purples = {[0.70, 0.50, 0.80]; ...
                     [0.55, 0.30, 0.65]; ...
                     [0.40, 0.15, 0.55]};

    % Sensitivity-analysis colours
    col_reset = [0.55, 0.30, 0.65];       % LCM purple
    col_decay = [0.65, 0.25, 0.50];       % warm-shifted purple
    col_logBF = [0.35, 0.45, 0.70];       % slate blue

    % Behavioural block colours (schematic only)
    col_pla = [0.30, 0.50, 0.90];         % placebo blue
    col_noc = [0.8431, 0.2549, 0.2549];   % nocebo red

    % Control correlation colour
    col_ctrl = [0.6, 0.6, 0.6];           % gray

    % Short labels
    shortLabels = {'Study 1', 'Study 2', 'Study 3'};

    %% ================================================================
    %%  FIGURE 1:  Model Comparison  (4 rows × 3 columns)
    %% ================================================================
    fig1 = figure('Color', 'w', 'Position', [50 50 1600 1100]);

    % --- Pre-compute consistent y-limits across studies for each row ---
    globalMax_dBIC  = 0;
    globalMin_logBF = 0;  globalMax_logBF = 0;
    globalMin_reset = 0;  globalMax_reset = 0;
    globalMin_decay = 0;  globalMax_decay = 0;

    for s = 1:nStudies
        R = allStudyResults{s};
        bestBIC = min([R.BIC_RW, R.BIC_BCC, R.BIC_LCM], [], 2);
        dAll = [R.BIC_RW - bestBIC; R.BIC_BCC - bestBIC; R.BIC_LCM - bestBIC];
        globalMax_dBIC = max(globalMax_dBIC, prctile(dAll, 99));

        globalMin_logBF = min(globalMin_logBF, prctile(R.logBF_LCM_vs_RW, 1));
        globalMax_logBF = max(globalMax_logBF, prctile(R.logBF_LCM_vs_RW, 99));

        diff_reset = R.BIC_LCM_reset - R.BIC_LCM;
        globalMin_reset = min(globalMin_reset, prctile(diff_reset, 1));
        globalMax_reset = max(globalMax_reset, prctile(diff_reset, 99));

        diff_decay = R.BIC_LCM_decay - R.BIC_LCM;
        globalMin_decay = min(globalMin_decay, prctile(diff_decay, 1));
        globalMax_decay = max(globalMax_decay, prctile(diff_decay, 99));
    end
    pad = 0.10;
    yl_dBIC  = [-5, globalMax_dBIC * (1 + pad)];
    yl_logBF = [globalMin_logBF - 1, globalMax_logBF * (1 + pad)];
    yl_reset = [globalMin_reset * (1 + pad), globalMax_reset * (1 + pad) + 5];
    yl_decay = [globalMin_decay * (1 + pad), globalMax_decay * (1 + pad) + 5];

    % ---- Row 1: BIC model comparison ----
    for s = 1:nStudies
        subplot(4, 3, s); hold on;
        R = allStudyResults{s};

        bestBIC  = min([R.BIC_RW, R.BIC_BCC, R.BIC_LCM], [], 2);
        dBIC_RW  = R.BIC_RW  - bestBIC;
        dBIC_BCC = R.BIC_BCC - bestBIC;
        dBIC_LCM = R.BIC_LCM - bestBIC;

        al_goodplot(dBIC_RW,  1, 0.35, col_RW,  'bilateral', [], std(dBIC_RW)/1000);
        al_goodplot(dBIC_BCC, 2, 0.35, col_BCC, 'bilateral', [], std(dBIC_BCC)/1000);
        al_goodplot(dBIC_LCM, 3, 0.35, col_LCM, 'bilateral', [], std(dBIC_LCM)/1000);

        % Percentage-best annotations
        yTop = yl_dBIC(2) * 0.92;
        text(1, yTop, sprintf('%.0f%%', R.pctBest_RW), ...
            'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', ...
            'Color', col_RW);
        text(2, yTop, sprintf('%.0f%%', R.pctBest_BCC), ...
            'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', ...
            'Color', col_BCC);
        text(3, yTop, sprintf('%.0f%%', R.pctBest_LCM), ...
            'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', ...
            'Color', col_LCM);

        set(gca, 'XTick', 1:3, 'XTickLabel', {'RW', 'BCC', 'LCM'}, 'FontSize', 12);
        ylim(yl_dBIC);
        if s == 1; ylabel('\DeltaBIC (from best model)', 'FontSize', 13); end
        title(studyLabels{s}, 'FontSize', 14);
        grid on; box off;

        if s == 1
            text(-0.15, 1.08, 'A', 'Units', 'normalized', ...
                'FontSize', 18, 'FontWeight', 'bold');
        end
    end

    % ---- Row 2: Log Bayes Factor (multiple causes vs single cause) ----
    for s = 1:nStudies
        subplot(4, 3, 3 + s); hold on;
        R = allStudyResults{s};

        al_goodplot(R.logBF_LCM_vs_RW, 1, 0.5, col_logBF, 'bilateral', [], ...
            std(R.logBF_LCM_vs_RW)/1000);

        % Threshold line at logBF = 3 (strong evidence)
        yline(3, '--', 'Color', [0.8, 0.3, 0.3], 'LineWidth', 1.2);
        yline(0, ':', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.8);

        % Annotations
        pctStrong = 100 * mean(R.logBF_LCM_vs_RW > 3);
        text(1, yl_logBF(2) * 0.90, ...
            sprintf('%.0f%% logBF > 3', pctStrong), ...
            'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', ...
            'Color', col_logBF);

        % Mean alpha annotation
        text(1, yl_logBF(1) + 0.5, ...
            sprintf('\\alpha = %.2f', mean(R.alpha)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0.4, 0.4, 0.4]);

        set(gca, 'XTick', 1, 'XTickLabel', {'LCM'}, 'FontSize', 12);
        ylim(yl_logBF);
        if s == 1; ylabel('log Bayes Factor', 'FontSize', 13); end
        grid on; box off;

        if s == 1
            text(-0.15, 1.08, 'B', 'Units', 'normalized', ...
                'FontSize', 18, 'FontWeight', 'bold');
        end
    end

    % ---- Row 3: Block-reset sensitivity ----
    for s = 1:nStudies
        subplot(4, 3, 6 + s); hold on;
        R = allStudyResults{s};

        diff_reset = R.BIC_LCM_reset - R.BIC_LCM;
        al_goodplot(diff_reset, 1, 0.5, col_reset, 'bilateral', [], ...
            std(diff_reset)/1000);

        yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

        nReset = sum(R.BIC_LCM_reset < R.BIC_LCM);
        pctReset = 100 * nReset / R.nSub;
        text(1, yl_reset(2) * 0.85, sprintf('%.0f%% prefer reset', pctReset), ...
            'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', ...
            'Color', col_reset);

        set(gca, 'XTick', 1, 'XTickLabel', {'LCM'}, 'FontSize', 12);
        ylim(yl_reset);
        if s == 1; ylabel('\DeltaBIC (reset − carry-over)', 'FontSize', 13); end
        grid on; box off;

        if s == 1
            text(-0.15, 1.08, 'C', 'Units', 'normalized', ...
                'FontSize', 18, 'FontWeight', 'bold');
        end
    end

    % ---- Row 4: Temporal decay sensitivity ----
    for s = 1:nStudies
        subplot(4, 3, 9 + s); hold on;
        R = allStudyResults{s};

        diff_decay = R.BIC_LCM_decay - R.BIC_LCM;
        al_goodplot(diff_decay, 1, 0.5, col_decay, 'bilateral', [], ...
            std(diff_decay)/1000);

        yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

        nDecay = sum(R.BIC_LCM_decay < R.BIC_LCM);
        pctDecay = 100 * nDecay / R.nSub;
        text(1, yl_decay(2) * 0.85, sprintf('%.0f%% prefer decay', pctDecay), ...
            'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', ...
            'Color', col_decay);

        set(gca, 'XTick', 1, 'XTickLabel', {'LCM'}, 'FontSize', 12);
        ylim(yl_decay);
        if s == 1; ylabel('\DeltaBIC (decay − standard)', 'FontSize', 13); end
        grid on; box off;

        if s == 1
            text(-0.15, 1.08, 'D', 'Units', 'normalized', ...
                'FontSize', 18, 'FontWeight', 'bold');
        end
    end

    sgtitle('Computational Modeling Results', 'FontSize', 18, 'FontWeight', 'bold');

    %% ================================================================
    %%  FIGURE 2:  Cause Assignment Mechanism  (PNAS double-column)
    %%  Layout:  A (schematic) | B (individual diffs)     — top
    %%           C (control correlations)                  — middle
    %%           D (learning vs novel) | E (entropy)       — bottom
    %% ================================================================
    % PNAS double-column: 17.8 cm ≈ 504 pt.  Aspect ~1.25:1 for 5-panel.
    fig2 = figure('Color', 'w', 'Position', [50 50 1500 1100]);

    % --- Shared settings for PNAS style ---
    fontAx    = 11;   % axis tick labels
    fontLabel = 12;   % axis titles
    fontTitle = 13;   % panel titles
    fontPanel = 16;   % panel letter (A, B, C …)
    fontLeg   = 11;   % legend entries

    markers     = {'o', 's', 'd'};
    markerSizes = [30, 35, 40];

    % ---- Panel A: Schematic (top-left) ----
    ax_schem = axes('Position', [0.04, 0.68, 0.42, 0.28]); hold on;
    draw_schematic(ax_schem, col_pla, col_noc);
    text(-0.05, 1.10, 'A', 'Units', 'normalized', ...
        'FontSize', fontPanel, 'FontWeight', 'bold');

    % ---- Panel B: Individual differences scatter (top-right) ----
    ax_corr = axes('Position', [0.56, 0.70, 0.40, 0.26]); hold on;

    legendEntries = cell(nStudies, 1);
    hScatter = gobjects(nStudies, 1);

    for s = 1:nStudies
        R = allStudyResults{s};
        valid = ~isnan(R.dir_soft_reset) & ~isnan(R.dir_obs);
        x = R.dir_soft_reset(valid);
        y = R.dir_obs(valid);

        hScatter(s) = scatter(x, y, markerSizes(s), study_purples{s}, ...
            markers{s}, 'filled', 'MarkerFaceAlpha', 0.55);

        % OLS regression line
        if numel(x) > 2
            p_coeff = polyfit(x, y, 1);
            x_range = linspace(min(x), max(x), 100);
            plot(x_range, polyval(p_coeff, x_range), '-', ...
                'Color', study_purples{s}, 'LineWidth', 1.5, ...
                'HandleVisibility', 'off');
        end

        stars = sigstars(R.pval_corr_soft);
        legendEntries{s} = sprintf('%s: \\rho = %.2f %s', ...
            shortLabels{s}, R.rho_soft, stars);
    end

    % Reference crosshairs
    xline(0, ':', 'Color', [0.7, 0.7, 0.7], 'HandleVisibility', 'off');
    yline(0, ':', 'Color', [0.7, 0.7, 0.7], 'HandleVisibility', 'off');

    xlabel('Model Cause-Assignment Shift', 'FontSize', fontLabel);
    ylabel('Observed VAS Diff (Noc − Pla)', 'FontSize', fontLabel);
    title('Individual Differences', 'FontSize', fontTitle);
    legend(hScatter, legendEntries, 'Location', 'northwest', 'FontSize', fontLeg, ...
        'Box', 'off');
    set(gca, 'FontSize', fontAx);
    grid on; box off;

    text(-0.05, 1.10, 'B', 'Units', 'normalized', ...
        'FontSize', fontPanel, 'FontWeight', 'bold');

    % ---- Panel C: Control correlations (middle row, full width) ----
    ax_ctrl = axes('Position', [0.08, 0.38, 0.88, 0.22]); hold on;

    % Check if control correlation fields exist
    hasCtrl = isfield(allStudyResults{1}, 'rho_eta');

    if hasCtrl
        predLabels = {'\eta_{RW}', '\kappa_{BCC}', '\alpha_{LCM}', 'LCM shift'};
        nPred = 4;

        rhoMat  = nan(nPred, nStudies);
        pvalMat = nan(nPred, nStudies);

        for s = 1:nStudies
            R = allStudyResults{s};
            rhoMat(1, s) = R.rho_eta;
            rhoMat(2, s) = R.rho_kappa;
            rhoMat(3, s) = R.rho_alpha;
            rhoMat(4, s) = R.rho_soft;

            pvalMat(1, s) = R.pval_eta;
            pvalMat(2, s) = R.pval_kappa;
            pvalMat(3, s) = R.pval_alpha;
            pvalMat(4, s) = R.pval_corr_soft;
        end

        barWidth = 0.2;
        studyOffsets = [-0.25, 0, 0.25];

        for p = 1:nPred
            for s = 1:nStudies
                xpos = p + studyOffsets(s);
                rho_val = rhoMat(p, s);
                if isnan(rho_val); continue; end

                if pvalMat(p, s) < .05
                    faceAlpha = 0.85;
                    edgeCol = study_purples{s} * 0.7;
                    edgeWidth = 1.2;
                else
                    faceAlpha = 0.25;
                    edgeCol = 'none';
                    edgeWidth = 0.5;
                end

                bar(xpos, rho_val, barWidth, 'FaceColor', study_purples{s}, ...
                    'EdgeColor', edgeCol, 'LineWidth', edgeWidth, ...
                    'FaceAlpha', faceAlpha);

                stars = sigstars(pvalMat(p, s));
                ypos = rho_val + sign(rho_val) * 0.03;
                if rho_val < 0; va = 'top'; else; va = 'bottom'; end
                text(xpos, ypos, stars, 'FontSize', 8, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', va);
            end
        end

        yline(0, '-', 'Color', [0.3, 0.3, 0.3], 'LineWidth', 0.8);
        xline(nPred - 0.5, ':', 'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.8);

        set(gca, 'XTick', 1:nPred, 'XTickLabel', predLabels, ...
            'FontSize', fontAx);
        ylabel('Spearman \rho', 'FontSize', fontLabel);
        title('Specificity: Control Correlations vs. LCM Cause-Assignment Shift', ...
            'FontSize', fontTitle);
        ylim([-0.5, 0.6]);

        % Compact legend (right side, inside axes)
        for s = 1:nStudies
            xleg = nPred + 0.50;
            yleg = 0.52 - (s - 1) * 0.08;
            patch([xleg-0.05, xleg+0.05, xleg+0.05, xleg-0.05], ...
                  [yleg-0.015, yleg-0.015, yleg+0.015, yleg+0.015], ...
                  study_purples{s}, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
            text(xleg + 0.10, yleg, shortLabels{s}, 'FontSize', fontLeg, ...
                'VerticalAlignment', 'middle');
        end

        text(0.99, 0.04, 'faded = n.s.', 'Units', 'normalized', ...
            'FontSize', 8, 'HorizontalAlignment', 'right', ...
            'Color', [0.5, 0.5, 0.5]);

        grid on; box off;

        text(-0.05, 1.10, 'C', 'Units', 'normalized', ...
            'FontSize', fontPanel, 'FontWeight', 'bold');
    else
        text(0.5, 0.5, 'Control correlations not computed (re-run run\_models)', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'FontSize', fontTitle, 'Color', [0.5, 0.5, 0.5]);
        axis off;
    end

    % ---- Panel D: Learning vs Novel group split (bottom-left) ----
    ax_group = axes('Position', [0.08, 0.06, 0.40, 0.24]); hold on;

    col_learn = [0.30, 0.70, 0.40];   % green
    col_novel = [0.85, 0.40, 0.30];   % red-orange

    if ~isempty(causeDiag)
        xPositions = [];
        groupData  = {};
        groupColor = {};
        xTickPos   = [];
        xTickLab   = {};
        idx = 0;

        for s = 1:nStudies
            D = causeDiag{s};
            behav = D.behav_effect;
            novel = D.assigned_novel;
            valid_mask = ~isnan(behav);

            learn_flag = valid_mask & ~novel;
            novel_flag = valid_mask & novel;

            idx = idx + 1;
            xPositions(end+1) = idx; %#ok<AGROW>
            d_learn = behav(learn_flag);
            if isempty(d_learn); d_learn = 0; end
            groupData{end+1}  = d_learn; %#ok<AGROW>
            groupColor{end+1} = col_learn; %#ok<AGROW>

            idx = idx + 1;
            xPositions(end+1) = idx; %#ok<AGROW>
            d_novel = behav(novel_flag);
            if isempty(d_novel); d_novel = 0; end
            groupData{end+1}  = d_novel; %#ok<AGROW>
            groupColor{end+1} = col_novel; %#ok<AGROW>

            xTickPos(end+1) = idx - 0.5; %#ok<AGROW>
            xTickLab{end+1} = shortLabels{s}; %#ok<AGROW>

            idx = idx + 0.5;
        end

        for g = 1:numel(groupData)
            d = groupData{g};
            if numel(d) > 2
                al_goodplot(d, xPositions(g), 0.35, groupColor{g}, 'bilateral', [], ...
                    std(d)/1000);
            elseif numel(d) >= 1
                bar(xPositions(g), mean(d), 0.35, 'FaceColor', groupColor{g}, ...
                    'EdgeColor', 'none', 'FaceAlpha', 0.6);
            end
        end

        yline(0, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
        set(gca, 'XTick', xTickPos, 'XTickLabel', xTickLab, 'FontSize', fontAx);
        ylabel('Observed VAS Diff (Noc − Pla)', 'FontSize', fontLabel);
        title('Behavioral Effect by Cause Assignment Type', 'FontSize', fontTitle);

        % Compact legend
        text(0.98, 0.95, 'Learning cause', 'Units', 'normalized', ...
            'FontSize', fontLeg, 'HorizontalAlignment', 'right', ...
            'Color', col_learn, 'FontWeight', 'bold');
        text(0.98, 0.87, 'Novel cause', 'Units', 'normalized', ...
            'FontSize', fontLeg, 'HorizontalAlignment', 'right', ...
            'Color', col_novel, 'FontWeight', 'bold');

        grid on; box off;
    else
        text(0.5, 0.5, 'Run diagnose\_causes first (pass as 3rd argument)', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'FontSize', fontTitle, 'Color', [0.5, 0.5, 0.5]);
        axis off;
    end

    text(-0.05, 1.10, 'D', 'Units', 'normalized', ...
        'FontSize', fontPanel, 'FontWeight', 'bold');

    % ---- Panel E: Entropy vs behavioral effect (bottom-right) ----
    ax_ent = axes('Position', [0.56, 0.06, 0.40, 0.24]); hold on;

    if ~isempty(causeDiag)
        hScatter_ent = gobjects(nStudies, 1);
        legendEntries_ent = cell(nStudies, 1);

        for s = 1:nStudies
            D = causeDiag{s};
            entropy_s = D.entropy_mix;
            behav_s   = D.behav_effect;

            valid_mask = ~isnan(entropy_s) & ~isnan(behav_s);
            x = entropy_s(valid_mask);
            y = behav_s(valid_mask);

            hScatter_ent(s) = scatter(x, y, markerSizes(s), study_purples{s}, ...
                markers{s}, 'filled', 'MarkerFaceAlpha', 0.55);

            if numel(x) > 2
                [rho_e, p_e] = corr(x, y, 'Type', 'Spearman');
                p_coeff = polyfit(x, y, 1);
                x_range = linspace(min(x), max(x), 100);
                plot(x_range, polyval(p_coeff, x_range), '-', ...
                    'Color', study_purples{s}, 'LineWidth', 1.5, ...
                    'HandleVisibility', 'off');
                legendEntries_ent{s} = sprintf('%s: \\rho = %.2f %s', ...
                    shortLabels{s}, rho_e, sigstars(p_e));
            else
                legendEntries_ent{s} = shortLabels{s};
            end
        end

        xlabel('Posterior Entropy on Mix Trials (bits)', 'FontSize', fontLabel);
        ylabel('Observed VAS Diff (Noc − Pla)', 'FontSize', fontLabel);
        title('Assignment Uncertainty vs. Behavioral Effect', 'FontSize', fontTitle);
        legend(hScatter_ent, legendEntries_ent, 'Location', 'northeast', ...
            'FontSize', fontLeg, 'Box', 'off');
        set(gca, 'FontSize', fontAx);
        grid on; box off;
    else
        text(0.5, 0.5, 'Run diagnose\_causes first (pass as 3rd argument)', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'FontSize', fontTitle, 'Color', [0.5, 0.5, 0.5]);
        axis off;
    end

    text(-0.05, 1.10, 'E', 'Units', 'normalized', ...
        'FontSize', fontPanel, 'FontWeight', 'bold');

end


%% ========================================================================
%  Helper: significance stars
%% ========================================================================
function s = sigstars(p)
    if isnan(p)
        s = '';
    elseif p < .001
        s = '***';
    elseif p < .01
        s = '**';
    elseif p < .05
        s = '*';
    else
        s = 'n.s.';
    end
end


%% ========================================================================
%  Helper: draw cause-assignment schematic  (Panel A of Figure 2)
%% ========================================================================
function draw_schematic(ax, col_pla, col_noc)
    axes(ax); %#ok<LAXES>

    xv = linspace(15, 85, 300);
    sigma = 6;

    % =========================
    %  TOP: Placebo block
    % =========================
    y_base_pla = 0.58;

    g_low_pla  = normpdf(xv, 40, sigma);
    g_low_pla  = g_low_pla / max(g_low_pla) * 0.28;

    g_high_pla = normpdf(xv, 70, sigma);
    g_high_pla = g_high_pla / max(g_high_pla) * 0.14;

    fill(xv, g_low_pla + y_base_pla, col_pla, ...
        'FaceAlpha', 0.45, 'EdgeColor', col_pla, 'LineWidth', 0.8);
    fill(xv, g_high_pla + y_base_pla, [0.65, 0.65, 0.65], ...
        'FaceAlpha', 0.30, 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.8);

    text(40, y_base_pla + 0.31, 'VAS 40', 'FontSize', 9, ...
        'HorizontalAlignment', 'center', 'Color', col_pla, 'FontWeight', 'bold');
    text(70, y_base_pla + 0.17, 'VAS 70', 'FontSize', 9, ...
        'HorizontalAlignment', 'center', 'Color', [0.45, 0.45, 0.45]);

    plot([50, 50], [y_base_pla - 0.02, y_base_pla + 0.35], '--k', 'LineWidth', 1.2);
    text(51, y_base_pla + 0.36, 'VAS 50', 'FontSize', 9, 'FontWeight', 'bold');

    annotation('textarrow', ...
        pos2norm(ax, [50, 43]), pos2norm_y(ax, [y_base_pla + 0.12, y_base_pla + 0.18]), ...
        'String', '', 'Color', col_pla, 'LineWidth', 1.5, 'HeadWidth', 8, 'HeadLength', 6);

    text(28, y_base_pla + 0.06, {'\rightarrow low-pain cause', '\rightarrow lower rating'}, ...
        'FontSize', 9, 'Color', col_pla, 'FontWeight', 'bold');

    text(15, y_base_pla + 0.38, 'Placebo block', 'FontSize', 12, ...
        'FontWeight', 'bold', 'Color', col_pla);

    % =========================
    %  BOTTOM: Nocebo block
    % =========================
    y_base_noc = 0.05;

    g_low_noc  = normpdf(xv, 30, sigma);
    g_low_noc  = g_low_noc / max(g_low_noc) * 0.14;

    g_high_noc = normpdf(xv, 60, sigma);
    g_high_noc = g_high_noc / max(g_high_noc) * 0.28;

    fill(xv, g_low_noc + y_base_noc, [0.65, 0.65, 0.65], ...
        'FaceAlpha', 0.30, 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.8);
    fill(xv, g_high_noc + y_base_noc, col_noc, ...
        'FaceAlpha', 0.45, 'EdgeColor', col_noc, 'LineWidth', 0.8);

    text(30, y_base_noc + 0.17, 'VAS 30', 'FontSize', 9, ...
        'HorizontalAlignment', 'center', 'Color', [0.45, 0.45, 0.45]);
    text(60, y_base_noc + 0.31, 'VAS 60', 'FontSize', 9, ...
        'HorizontalAlignment', 'center', 'Color', col_noc, 'FontWeight', 'bold');

    plot([50, 50], [y_base_noc - 0.02, y_base_noc + 0.35], '--k', 'LineWidth', 1.2);
    text(51, y_base_noc + 0.36, 'VAS 50', 'FontSize', 9, 'FontWeight', 'bold');

    annotation('textarrow', ...
        pos2norm(ax, [50, 57]), pos2norm_y(ax, [y_base_noc + 0.12, y_base_noc + 0.18]), ...
        'String', '', 'Color', col_noc, 'LineWidth', 1.5, 'HeadWidth', 8, 'HeadLength', 6);

    text(62, y_base_noc + 0.06, {'\rightarrow high-pain cause', '\rightarrow higher rating'}, ...
        'FontSize', 9, 'Color', col_noc, 'FontWeight', 'bold');

    text(15, y_base_noc + 0.38, 'Nocebo block', 'FontSize', 12, ...
        'FontWeight', 'bold', 'Color', col_noc);

    % --- Shared VAS axis at bottom ---
    plot([20, 80], [y_base_noc - 0.04, y_base_noc - 0.04], 'k-', 'LineWidth', 1);
    for tick = [20, 30, 40, 50, 60, 70, 80]
        plot([tick, tick], [y_base_noc - 0.06, y_base_noc - 0.04], 'k-', 'LineWidth', 0.8);
        text(tick, y_base_noc - 0.09, num2str(tick), 'FontSize', 8, ...
            'HorizontalAlignment', 'center');
    end
    text(50, y_base_noc - 0.14, 'Pain intensity (VAS)', 'FontSize', 10, ...
        'HorizontalAlignment', 'center');

    plot([18, 82], [0.53, 0.53], ':', 'Color', [0.75, 0.75, 0.75], 'LineWidth', 0.8);

    xlim([12, 88]);
    ylim([-0.15, 1.02]);
    axis off;
    title('Cause Assignment Logic', 'FontSize', 13);
end


%% ========================================================================
%  Helper: convert axis x-coords to normalised figure coords (for arrows)
%% ========================================================================
function xn = pos2norm(ax, xdata)
    xl = ax.XLim;
    pos = ax.Position;
    xn = pos(1) + (xdata - xl(1)) / (xl(2) - xl(1)) * pos(3);
end

function yn = pos2norm_y(ax, ydata)
    yl = ax.YLim;
    pos = ax.Position;
    yn = pos(2) + (ydata - yl(1)) / (yl(2) - yl(1)) * pos(4);
end
