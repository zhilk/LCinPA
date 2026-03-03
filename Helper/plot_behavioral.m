function plot_behavioral(T, studyLabel)
%PLOT_BEHAVIORAL  Behavioral VAS rating triptych for one LCPA study.
%
%  INPUTS:
%    T           – table for one study (from readtable on the CSV)
%    studyLabel  – string label for sgtitle (e.g. 'Study 1 (Behavioural)')
%
%  Produces one figure with three subplots:
%    Left   – Placebo block rainclouds  (low / mix / high)
%    Centre – Mix trial violins with slope-coded subject lines
%    Right  – Nocebo block rainclouds   (low / mix / high)
%
%  Requires: raincloud_plot.m, al_goodplot.m on the MATLAB path.

%% ---- Compute subject-level means for test phase ----

test = T(strcmp(T.Phase, 'test'), :);

XPC  = subjectMeans(test, 'placebo', 'low');     % placebo low
X1PC = subjectMeans(test, 'placebo', 'mix');     % placebo mix
X2PC = subjectMeans(test, 'placebo', 'high');    % placebo high

XNC  = subjectMeans(test, 'nocebo',  'low');     % nocebo low
X1NC = subjectMeans(test, 'nocebo',  'mix');     % nocebo mix
X2NC = subjectMeans(test, 'nocebo',  'high');    % nocebo high

nSub = numel(X1PC);

%% ---- Adaptive y-limits for raincloud panels ----
if nSub <= 50
    ylimits = [-0.02, 0.03];
else
    ylimits = [-0.025, 0.04];
end

%% ---- Colour palette ----
soft_blue1 = [0.7, 0.85, 1];           % light blue  (low)
soft_blue2 = [0.3, 0.5, 0.9];          % medium blue (mix)
soft_blue3 = [0.2, 0.3, 0.7];          % dark blue   (high)

soft_red1  = [0.89, 0.55, 0.55];       % light red   (low)
soft_red2  = [0.8431, 0.2549, 0.2549]; % medium red  (mix)
soft_red3  = [0.5, 0, 0];              % dark red    (high)

slope_pos = soft_blue2;                 % nocebo > placebo
slope_neg = [0.75, 0.3, 0.3];          % placebo > nocebo

% Shared raincloud options:
%   box_col_match = 0  → keeps bxcl black for edges, median, whiskers
%   bxfacecl      = color → fills the box with the cloud color
%   bxfacealpha   = 0.5 → semi-transparent box fill
%   cloud_edge_col= 'none' → no outline on the density cloud
rc_alpha = 0.6;     % cloud face alpha
bx_alpha = 0.8;     % box face alpha

%% ======================================================================
%  Single triptych figure
%  ======================================================================
figure('Color','w','Position',[100 100 1400 500]);

% --- Subplot 1: Placebo Block rainclouds ---
subplot(1,3,1);
raincloud_plot(XPC, 'box_on',1,'color',soft_blue1,'alpha',rc_alpha, ...
    'box_dodge',1,'box_dodge_amount',.15,'dot_dodge_amount',.15, ...
    'box_col_match',0,'bxfacecl',soft_blue1,'bxfacealpha',bx_alpha, ...
    'cloud_edge_col','none','line_width',0.7);
hold on;
raincloud_plot(X1PC,'box_on',1,'color',soft_blue2,'alpha',rc_alpha, ...
    'box_dodge',1,'box_dodge_amount',.35,'dot_dodge_amount',.35, ...
    'box_col_match',0,'bxfacecl',soft_blue2,'bxfacealpha',bx_alpha, ...
    'cloud_edge_col','none','line_width',0.7);
raincloud_plot(X2PC,'box_on',1,'color',soft_blue3,'alpha',rc_alpha, ...
    'box_dodge',1,'box_dodge_amount',.55,'dot_dodge_amount',.55, ...
    'box_col_match',0,'bxfacecl',soft_blue3,'bxfacealpha',bx_alpha, ...
    'cloud_edge_col','none','line_width',0.7);
set(gca, 'XLim',[-20 120], 'YLim',ylimits, 'GridAlpha',0.3);
xlabel('VAS Ratings','FontSize',14);
title('Placebo Block','FontSize',14);
grid on; box off;

% --- Subplot 2: Mix trial violins + slope lines ---
subplot(1,3,2);
hold on;

placebo_x = 0.5;
nocebo_x  = 0.55;

al_goodplot(X1PC, placebo_x, 0.3, soft_blue2, 'left',  [], std(X1PC)/1000);
al_goodplot(X1NC, nocebo_x,  0.3, soft_red2,  'right', [], std(X1NC)/1000);

for i = 1:nSub
    if X1NC(i) > X1PC(i)
        lc = slope_pos;
    else
        lc = slope_neg;
    end
    plot([placebo_x nocebo_x], [X1PC(i) X1NC(i)], '-', ...
        'Color', lc, 'LineWidth', 0.75);
end

hold off;
ylim([0 100]);
xlim([0.2 0.9]);
set(gca, 'XTickLabel', []);
ylabel('VAS Ratings','FontSize',14);
title('Mix Trials','FontSize',14);
grid on; box off;

% --- Subplot 3: Nocebo Block rainclouds ---
subplot(1,3,3);
raincloud_plot(XNC, 'box_on',1,'color',soft_red1,'alpha',rc_alpha, ...
    'box_dodge',1,'box_dodge_amount',.15,'dot_dodge_amount',.15, ...
    'box_col_match',0,'bxfacecl',soft_red1,'bxfacealpha',bx_alpha, ...
    'cloud_edge_col','none','line_width',0.7);
hold on;
raincloud_plot(X1NC,'box_on',1,'color',soft_red2,'alpha',rc_alpha, ...
    'box_dodge',1,'box_dodge_amount',.35,'dot_dodge_amount',.35, ...
    'box_col_match',0,'bxfacecl',soft_red2,'bxfacealpha',bx_alpha, ...
    'cloud_edge_col','none','line_width',0.7);
raincloud_plot(X2NC,'box_on',1,'color',soft_red3,'alpha',rc_alpha, ...
    'box_dodge',1,'box_dodge_amount',.55,'dot_dodge_amount',.55, ...
    'box_col_match',0,'bxfacecl',soft_red3,'bxfacealpha',bx_alpha, ...
    'cloud_edge_col','none','line_width',0.7);
set(gca, 'XLim',[-20 120], 'YLim',ylimits, 'GridAlpha',0.3);
xlabel('VAS Ratings','FontSize',14);
title('Nocebo Block','FontSize',14);
grid on; box off;

sgtitle(studyLabel, 'FontSize', 18, 'FontWeight', 'bold');

end


%% ========================================================================
%  LOCAL HELPER
%% ========================================================================
function m = subjectMeans(test, blockName, trialType)
%SUBJECTMEANS  Per-subject mean VASRating for a Block x trial type.
%
%  trialType: 'low'/'high' filter on CueAssociation; 'mix' on VisualCategory.

    rows = strcmp(test.Block, blockName);

    if strcmp(trialType, 'mix')
        rows = rows & strcmp(test.VisualCategory, 'mix');
    else
        rows = rows & strcmp(test.CueAssociation, trialType);
    end

    sub = test(rows, :);

    [G, subIDs] = findgroups(sub.SubID);
    m = splitapply(@(x) mean(x, 'omitnan'), sub.VASRating, G);

    [~, ord] = sort(subIDs);
    m = m(ord);
end
