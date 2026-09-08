function plot_conditioning(behavioralData, expIdx, opts)

%  (1) per-subject High vs Low means, (2) time course during conditioning.
%  ugly at the moment - can be fixed if to be included. 

    arguments
        behavioralData cell
        expIdx (1,1) double = 1
        opts.ExcludeNoResponse (1,1) logical = true
        opts.ExcludeCatch      (1,1) logical = true
    end

    T = behavioralData{expIdx};
    keep = strcmp(T.Phase,'conditioning');
    if opts.ExcludeNoResponse && ismember('VASResponse', T.Properties.VariableNames)
        keep = keep & T.VASResponse == 1;
    end
    if opts.ExcludeCatch && ismember('CatchTrial', T.Properties.VariableNames)
        keep = keep & T.CatchTrial == 0;
    end
    C = sortrows(T(keep,:), {'SubID','TrialInPhase'});   % order for the time course

    subs = unique(C.SubID);
    nS   = numel(subs);

    % ---------- per-subject means ----------
    hMean = nan(nS,1); lMean = nan(nS,1);
    for s = 1:nS
        r = C.SubID == subs(s);
        hMean(s) = mean(C.VASRating(r & strcmp(C.CueAssociation,'high')), 'omitnan');
        lMean(s) = mean(C.VASRating(r & strcmp(C.CueAssociation,'low')),  'omitnan');
    end

    figure('Color','w','Name',sprintf('Exp %d: per-subject means', expIdx)); hold on
    for s = 1:nS
        plot([1 2], [hMean(s) lMean(s)], '-', 'Color', [.7 .7 .7 .4]);
    end
    plot(1, hMean, 'o', 'MarkerFaceColor', [.2 .4 .8], 'MarkerEdgeColor','none');
    plot(2, lMean, 'o', 'MarkerFaceColor', [.85 .3 .3], 'MarkerEdgeColor','none');
    errorbar([1 2], [mean(hMean,'omitnan') mean(lMean,'omitnan')], ...
             [std(hMean,'omitnan') std(lMean,'omitnan')] ./ sqrt(nS), ...
             'k-', 'LineWidth', 2, 'CapSize', 12, 'Marker','s', 'MarkerFaceColor','k');
    xlim([.5 2.5]); xticks([1 2]); xticklabels({'High','Low'});
    ylabel('Mean VAS rating (conditioning)'); box off
    title(sprintf('Exp %d: per-subject High vs Low (n=%d)', expIdx, nS));

    % ---------- time course (ordinal within category) ----------
    maxH = max(arrayfun(@(s) sum(C.SubID==s & strcmp(C.CueAssociation,'high')), subs));
    maxF = max(arrayfun(@(s) sum(C.SubID==s & strcmp(C.CueAssociation,'low')),  subs));
    Hmat = nan(nS, maxH); Lmat = nan(nS, maxF);
    for s = 1:nS
        rs = C(C.SubID==subs(s), :);                       % already time-sorted
        hv = rs.VASRating(strcmp(rs.CueAssociation,'high'));
        lv = rs.VASRating(strcmp(rs.CueAssociation,'low'));
        Hmat(s,1:numel(hv)) = hv;
        Lmat(s,1:numel(lv)) = lv;
    end
    mH = mean(Hmat,1,'omitnan'); seH = std(Hmat,0,1,'omitnan')./sqrt(sum(~isnan(Hmat),1));
    mL = mean(Lmat,1,'omitnan'); seL = std(Lmat,0,1,'omitnan')./sqrt(sum(~isnan(Lmat),1));

    figure('Color','w','Name',sprintf('Exp %d: conditioning time course', expIdx)); hold on
    errorbar(1:maxH, mH, seH, '-o', 'Color',[.2 .4 .8], 'MarkerFaceColor',[.2 .4 .8], ...
             'CapSize',0, 'DisplayName','High');
    errorbar(1:maxF, mL, seL, '-o', 'Color',[.85 .3 .3], 'MarkerFaceColor',[.85 .3 .3], ...
             'CapSize',0, 'DisplayName','Low');
    xlabel('Trial number within category (conditioning)');
    ylabel('Mean VAS rating \pm SEM'); box off
    legend('Location','best');
    title(sprintf('Exp %d: conditioning time course', expIdx));
end