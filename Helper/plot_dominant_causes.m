function plot_dominant_causes(domCause, studyLabels)
%PLOT_DOMINANT_CAUSES  Visualize dominant cause structure from block-reset LCM.
%
%   plot_dominant_causes(domCause, studyLabels)
%
%   Produces one figure with two panels:
%
%   Panel A — "Cause Portrait": paired dot plot of dominant cause US means.
%   Panel B — "Posterior Weight Asymmetry": grouped bar chart.
%             Bar fill  = block mix-trial colour (matching behavioral figs).
%             Bar edge  = cause-affiliated conditioning colour (thick).
%
%   INPUTS:
%     domCause    – {nStudies×1} cell of structs from analyze_dominant_causes()
%     studyLabels – {nStudies×1} cell of study name strings
%
%   Requires: al_goodplot.m on the MATLAB path.

nStudies = numel(domCause);

% ---- Colour palette (matching behavioral al_goodplot figure exactly) ----
% al_goodplot uses 0.97*col for patch fill, facealpha 0.8 for box body
% Fill = 0.97 * soft_blue2 / soft_red2  (matching al_goodplot box body)
col_pla      = 0.97 * [0.30, 0.50, 0.90];          % placebo mix (al_goodplot)
col_noc      = 0.97 * [0.8431, 0.2549, 0.2549];    % nocebo mix  (al_goodplot)
barAlpha     = 0.8;                                  % matches al_goodplot box facealpha

% Edge = cause-affiliated conditioning colours (soft_blue1/3, soft_red1/3)
col_low_pla  = [0.70, 0.85, 1.00];          % placebo low-pain cause (light blue)
col_high_pla = [0.20, 0.30, 0.70];          % placebo high-pain cause (dark blue)
col_low_noc  = [0.89, 0.55, 0.55];          % nocebo low-pain cause  (light red)
col_high_noc = [0.50, 0.00, 0.00];          % nocebo high-pain cause (dark red)

col_other    = [0.75, 0.75, 0.75];          % other causes (gray)

% Edge line width for cause affiliation (thick!)
edgeLW       = 4.5;

% Short labels with study description
shortLabels = {'Study 1 (Behavioral)', 'Study 2 (Control)', 'Study 3 (fMRI)'};
if numel(studyLabels) < 3
    shortLabels = shortLabels(1:numel(studyLabels));
end

% Font sizes
fontAx    = 11;
fontLabel = 13;
fontTitle = 14;
fontPanel = 18;

% ================================================================
%  FIGURE: Dominant Cause Structure  (2 rows x nStudies columns)
% ================================================================
fig = figure('Color', 'w', 'Position', [80 80 500*nStudies 900]);

% ────────────────────────────────────────────────────────────────
%  ROW 1 — Panel A: Cause Portrait (paired dot plot)
% ────────────────────────────────────────────────────────────────
for s = 1:nStudies
    ax = subplot(2, nStudies, s); hold on;

    D = domCause{s};
    nSub = size(D.top2_us_pla, 1);

    % Sort subjects by their mean placebo cause US for visual clarity
    mean_pla = mean(D.top2_us_pla, 2, 'omitnan');
    [~, sortOrder] = sort(mean_pla);

    % -- Placebo block (left side: x = 1 and 2) --
    us_pla = D.top2_us_pla(sortOrder, :);
    % -- Nocebo block (right side: x = 4 and 5) --
    us_noc = D.top2_us_noc(sortOrder, :);

    % Jitter for subjects
    jit_pla = linspace(-0.15, 0.15, nSub);
    jit_noc = linspace(-0.15, 0.15, nSub);

    % Plot individual subject lines connecting their two dominant causes
    for i = 1:nSub
        % Placebo pair (blue connector)
        if ~any(isnan(us_pla(i,:)))
            plot([1 + jit_pla(i), 2 + jit_pla(i)], us_pla(i,:), '-', ...
                'Color', [col_pla, 0.12], 'LineWidth', 0.8);
        end
        % Nocebo pair (red connector)
        if ~any(isnan(us_noc(i,:)))
            plot([4 + jit_noc(i), 5 + jit_noc(i)], us_noc(i,:), '-', ...
                'Color', [col_noc, 0.12], 'LineWidth', 0.8);
        end
    end

    % Scatter: cause colours (light = low, dark = high)
    for i = 1:nSub
        if ~isnan(us_pla(i,1))
            scatter(1 + jit_pla(i), us_pla(i,1), 18, col_low_pla, 'o', ...
                'filled', 'MarkerFaceAlpha', 0.45);
        end
        if ~isnan(us_pla(i,2))
            scatter(2 + jit_pla(i), us_pla(i,2), 18, col_high_pla, 'o', ...
                'filled', 'MarkerFaceAlpha', 0.45);
        end
        if ~isnan(us_noc(i,1))
            scatter(4 + jit_noc(i), us_noc(i,1), 18, col_low_noc, 'o', ...
                'filled', 'MarkerFaceAlpha', 0.45);
        end
        if ~isnan(us_noc(i,2))
            scatter(5 + jit_noc(i), us_noc(i,2), 18, col_high_noc, 'o', ...
                'filled', 'MarkerFaceAlpha', 0.45);
        end
    end

    % Group means with error bars (SE)
    means = [mean(us_pla(:,1),'omitnan'), mean(us_pla(:,2),'omitnan'), ...
             NaN, ...
             mean(us_noc(:,1),'omitnan'), mean(us_noc(:,2),'omitnan')];
    sems  = [std(us_pla(:,1),'omitnan')/sqrt(sum(~isnan(us_pla(:,1)))), ...
             std(us_pla(:,2),'omitnan')/sqrt(sum(~isnan(us_pla(:,2)))), ...
             NaN, ...
             std(us_noc(:,1),'omitnan')/sqrt(sum(~isnan(us_noc(:,1)))), ...
             std(us_noc(:,2),'omitnan')/sqrt(sum(~isnan(us_noc(:,2))))];

    xpos = [1, 2, 3, 4, 5];
    for j = [1, 2, 4, 5]
        if ~isnan(means(j))
            errorbar(xpos(j), means(j), sems(j), 'o', ...
                'MarkerSize', 10, 'LineWidth', 2.0, ...
                'Color', 'k', 'MarkerFaceColor', 'k', ...
                'CapSize', 8);
        end
    end

    % Connect means with thick lines
    plot([1, 2], means([1, 2]), '-k', 'LineWidth', 2.0);
    plot([4, 5], means([4, 5]), '-k', 'LineWidth', 2.0);

    % VAS 50 reference line (dashed)
    yline(50, '--', 'Color', [0.3, 0.3, 0.3], 'LineWidth', 1.2);
    text(5.7, 50, 'VAS 50', 'FontSize', 9, 'Color', [0.3, 0.3, 0.3], ...
        'VerticalAlignment', 'bottom');

    % Separator between blocks
    xline(3, '-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);

    % Block labels above the data
    text(1.5, 83, 'Placebo', 'FontSize', 12, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'Color', col_pla);
    text(4.5, 83, 'Nocebo', 'FontSize', 12, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'Color', col_noc);

    % X-axis
    set(gca, 'XTick', [1, 2, 4, 5], ...
        'XTickLabel', {'Low', 'High', 'Low', 'High'}, ...
        'FontSize', fontAx);
    xlabel('Cause Affiliation', 'FontSize', fontLabel);
    xlim([0.3, 6.2]);
    ylim([15, 87]);
    if s == 1
        ylabel('Inferred US Mean (VAS)', 'FontSize', fontLabel);
    end
    title(shortLabels{s}, 'FontSize', fontTitle);
    grid on; box off;

    if s == 1
        text(-0.12, 1.06, 'A', 'Units', 'normalized', ...
            'FontSize', fontPanel, 'FontWeight', 'bold');
    end
end

% ────────────────────────────────────────────────────────────────
%  ROW 2 — Panel B: Posterior Weight Asymmetry (grouped bars)
% ────────────────────────────────────────────────────────────────
for s = 1:nStudies
    ax = subplot(2, nStudies, nStudies + s); hold on;

    D = domCause{s};

    % Bar positions
    %   x=1: Pla-Low,  x=2: Pla-High
    %   x=4: Noc-Low,  x=5: Noc-High
    barW = 0.6;
    halfW = barW / 2;

    % Mean values
    m_pla_low   = mean(D.postw_low_pla,  'omitnan');
    m_pla_high  = mean(D.postw_high_pla, 'omitnan');
    m_noc_low   = mean(D.postw_low_noc,  'omitnan');
    m_noc_high  = mean(D.postw_high_noc, 'omitnan');
    m_pla_other = mean(D.postw_other_pla, 'omitnan');
    m_noc_other = mean(D.postw_other_noc, 'omitnan');

    % SE values
    n_pla = sum(~isnan(D.postw_low_pla));
    n_noc = sum(~isnan(D.postw_low_noc));
    se_pla_low   = std(D.postw_low_pla,  'omitnan') / sqrt(n_pla);
    se_pla_high  = std(D.postw_high_pla, 'omitnan') / sqrt(n_pla);
    se_noc_low   = std(D.postw_low_noc,  'omitnan') / sqrt(n_noc);
    se_noc_high  = std(D.postw_high_noc, 'omitnan') / sqrt(n_noc);

    % Draw bars: mix-trial fill + thick cause-affiliated edge
    draw_bar_edge(1, m_pla_low,  halfW, col_pla, barAlpha, col_low_pla,  edgeLW);
    draw_bar_edge(2, m_pla_high, halfW, col_pla, barAlpha, col_high_pla, edgeLW);
    draw_bar_edge(4, m_noc_low,  halfW, col_noc, barAlpha, col_low_noc,  edgeLW);
    draw_bar_edge(5, m_noc_high, halfW, col_noc, barAlpha, col_high_noc, edgeLW);

    % Gray "other" stacked on top of each bar
    draw_bar_stacked(1, m_pla_low,  m_pla_other, halfW, col_other, 0.35);
    draw_bar_stacked(2, m_pla_high, m_pla_other, halfW, col_other, 0.35);
    draw_bar_stacked(4, m_noc_low,  m_noc_other, halfW, col_other, 0.35);
    draw_bar_stacked(5, m_noc_high, m_noc_other, halfW, col_other, 0.35);

    % Error bars (on main bars only)
    errorbar(1, m_pla_low,  se_pla_low,  'k', 'LineWidth', 1.5, 'CapSize', 6);
    errorbar(2, m_pla_high, se_pla_high, 'k', 'LineWidth', 1.5, 'CapSize', 6);
    errorbar(4, m_noc_low,  se_noc_low,  'k', 'LineWidth', 1.5, 'CapSize', 6);
    errorbar(5, m_noc_high, se_noc_high, 'k', 'LineWidth', 1.5, 'CapSize', 6);

    % Individual data points: dark grey for visibility on any bar fill
    jit_w = 0.22;
    dot_col = [0.25, 0.25, 0.25];   % uniform dark grey
    dot_alpha = 0.30;
    valid_pla = ~isnan(D.postw_low_pla);
    valid_noc = ~isnan(D.postw_low_noc);

    rng(42 + s);  % reproducible jitter per study
    jit1 = (rand(sum(valid_pla), 1) - 0.5) * jit_w;
    scatter(1 + jit1, D.postw_low_pla(valid_pla),  10, dot_col, ...
        'filled', 'MarkerFaceAlpha', dot_alpha);
    jit2 = (rand(sum(valid_pla), 1) - 0.5) * jit_w;
    scatter(2 + jit2, D.postw_high_pla(valid_pla), 10, dot_col, ...
        'filled', 'MarkerFaceAlpha', dot_alpha);
    jit3 = (rand(sum(valid_noc), 1) - 0.5) * jit_w;
    scatter(4 + jit3, D.postw_low_noc(valid_noc),  10, dot_col, ...
        'filled', 'MarkerFaceAlpha', dot_alpha);
    jit4 = (rand(sum(valid_noc), 1) - 0.5) * jit_w;
    scatter(5 + jit4, D.postw_high_noc(valid_noc), 10, dot_col, ...
        'filled', 'MarkerFaceAlpha', dot_alpha);

    % Significance brackets
    ymax_bars = max([m_pla_low + m_pla_other + se_pla_low, ...
                     m_pla_high + m_pla_other + se_pla_high, ...
                     m_noc_low + m_noc_other + se_noc_low, ...
                     m_noc_high + m_noc_other + se_noc_high]);
    bracket_y = ymax_bars + 0.05;

    % Placebo bracket (low > high?)
    shift_pla = D.postw_low_pla - D.postw_high_pla;
    valid_sp = ~isnan(shift_pla);
    if sum(valid_sp) >= 5
        [p_sp] = signrank(shift_pla(valid_sp));
        stars_pla = sigstars(p_sp);
        plot([1, 1, 2, 2], [bracket_y - 0.01, bracket_y, bracket_y, bracket_y - 0.01], ...
            'k-', 'LineWidth', 1.2);
        text(1.5, bracket_y + 0.015, stars_pla, 'FontSize', 12, ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end

    % Nocebo bracket (high > low?)
    shift_noc = D.postw_high_noc - D.postw_low_noc;
    valid_sn = ~isnan(shift_noc);
    if sum(valid_sn) >= 5
        [p_sn] = signrank(shift_noc(valid_sn));
        stars_noc = sigstars(p_sn);
        plot([4, 4, 5, 5], [bracket_y - 0.01, bracket_y, bracket_y, bracket_y - 0.01], ...
            'k-', 'LineWidth', 1.2);
        text(4.5, bracket_y + 0.015, stars_noc, 'FontSize', 12, ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end

    % Block separator
    xline(3, '-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);

    % Block labels above bars
    text(1.5, bracket_y + 0.07, 'Placebo', 'FontSize', 12, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'Color', col_pla);
    text(4.5, bracket_y + 0.07, 'Nocebo', 'FontSize', 12, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'Color', col_noc);

    % Axis formatting
    set(gca, 'XTick', [1, 2, 4, 5], ...
        'XTickLabel', {'Low', 'High', 'Low', 'High'}, ...
        'FontSize', fontAx);
    xlabel('Cause Affiliation', 'FontSize', fontLabel);
    xlim([0.3, 6.2]);
    ylim([0, bracket_y + 0.12]);
    if s == 1
        ylabel('Posterior Weight on Mix Trials', 'FontSize', fontLabel);
    end
    title(shortLabels{s}, 'FontSize', fontTitle);
    grid on; box off;

    % Panel label (only on first subplot)
    if s == 1
        text(-0.12, 1.06, 'B', 'Units', 'normalized', ...
            'FontSize', fontPanel, 'FontWeight', 'bold');
    end
end

sgtitle('Dominant Cause Structure', 'FontSize', 18, 'FontWeight', 'bold');

% ── Minimal legend centred at bottom ──
legend_y  = 0.02;
legend_fs = 10;
sw = 0.015; sh = 0.012;

annotation('rectangle', [0.48, legend_y + 0.003, sw, sh], ...
    'FaceColor', col_other, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
annotation('textbox', [0.48 + sw + 0.005, legend_y - 0.002, 0.12, 0.025], ...
    'String', 'Other causes', 'FontSize', legend_fs, ...
    'EdgeColor', 'none', 'VerticalAlignment', 'middle', ...
    'Color', [0.5, 0.5, 0.5]);

end


%% ========================================================================
%  Helper: draw a bar using patch (no BaseValue artifact)
%% ========================================================================
function draw_bar(xc, height, halfW, col, faceAlpha)
    if isnan(height) || height == 0; return; end
    patch([xc-halfW, xc+halfW, xc+halfW, xc-halfW], ...
          [0, 0, height, height], col, ...
          'FaceAlpha', faceAlpha, 'EdgeColor', 'none');
end


%% ========================================================================
%  Helper: draw a bar with a thick coloured edge
%% ========================================================================
function draw_bar_edge(xc, height, halfW, fillCol, faceAlpha, edgeCol, edgeLW)
    if isnan(height) || height == 0; return; end
    patch([xc-halfW, xc+halfW, xc+halfW, xc-halfW], ...
          [0, 0, height, height], fillCol, ...
          'FaceAlpha', faceAlpha, 'EdgeColor', edgeCol, 'LineWidth', edgeLW);
end


%% ========================================================================
%  Helper: draw a stacked segment on top of a base bar using patch
%% ========================================================================
function draw_bar_stacked(xc, base_height, seg_height, halfW, col, faceAlpha)
    if isnan(seg_height) || seg_height == 0; return; end
    if isnan(base_height); base_height = 0; end
    patch([xc-halfW, xc+halfW, xc+halfW, xc-halfW], ...
          [base_height, base_height, base_height + seg_height, base_height + seg_height], ...
          col, 'FaceAlpha', faceAlpha, 'EdgeColor', 'none');
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
