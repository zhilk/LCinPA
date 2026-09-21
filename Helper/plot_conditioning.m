%% conditioning: overlay experiments, faceted by BLOCK TYPE (nocebo/placebo)
expList           = 1:3;
ExcludeNoResponse = true;
ExcludeCatch      = true;
runVar            = 'Block';                         % <-- manipulation, not BlockOrder
expCols           = [.20 .40 .80; .85 .30 .30; .30 .65 .35];
cues              = {'high','low'};

% ---------- long table across experiments ----------
LT = table();
for e = expList
    T    = behavioralData{e};
    keep = strcmp(T.Phase,'conditioning');
    if ExcludeNoResponse && ismember('VASResponse',T.Properties.VariableNames)
        keep = keep & T.VASResponse == 1;
    end
    if ExcludeCatch && ismember('CatchTrial',T.Properties.VariableNames)
        keep = keep & T.CatchTrial == 0;
    end
    Te = sortrows(T(keep,:), {'SubID','Block','TrialInBlock'});
    Te.exp   = repmat(e, height(Te), 1);
    Te.runid = string(Te.(runVar));                  % block type as string

    % ordinal trial within subject x block x cue (resets per block)
    g   = findgroups(Te.SubID, Te.runid, Te.CueAssociation);
    ord = zeros(height(Te),1);
    for k = 1:max(g)
        idx = find(g==k);  ord(idx) = 1:numel(idx);
    end
    Te.ord = ord;

    LT = [LT; Te(:, {'exp','runid','SubID','CueAssociation','VASRating','ord'})];
end

runs = unique(LT.runid);        % ["nocebo","placebo"]
nRun = numel(runs);

% ---------- per-subject / group means ----------
figure('Color','w','Name','Conditioning: High vs Low means');
for ri = 1:nRun
    subplot(1,nRun,ri); hold on
    for ei = 1:numel(expList)
        e = expList(ei);
        for ci = 1:2
            sel  = LT.exp==e & LT.runid==runs(ri) & strcmp(LT.CueAssociation,cues{ci});
            subs = unique(LT.SubID(sel));
            if isempty(subs), continue; end
            pm   = arrayfun(@(s) mean(LT.VASRating(sel & LT.SubID==s),'omitnan'), subs);
            xpos = ci + (ei-2)*0.12;
            errorbar(xpos, mean(pm,'omitnan'), std(pm,'omitnan')/sqrt(numel(subs)), ...
                     'o', 'Color', expCols(ei,:), 'MarkerFaceColor', expCols(ei,:), ...
                     'CapSize', 8, 'LineWidth', 1.2, 'HandleVisibility','off');
        end
    end
    xlim([.5 2.5]); xticks([1 2]); xticklabels({'High','Low'});
    ylabel('Mean VAS'); box off; title(runs(ri));
    if ri==nRun
        for ei = 1:numel(expList)
            plot(nan,nan,'o','Color',expCols(ei,:),'MarkerFaceColor',expCols(ei,:), ...
                 'DisplayName',sprintf('Exp %d',expList(ei)));
        end
        legend('show','Location','best');
    end
end

% ---------- time course (ordinal within category, reset per block) ----------
figure('Color','w','Name','Conditioning: time course');
for ri = 1:nRun
    subplot(1,nRun,ri); hold on
    for ei = 1:numel(expList)
        e = expList(ei);
        for ci = 1:2
            sel  = LT.exp==e & LT.runid==runs(ri) & strcmp(LT.CueAssociation,cues{ci});
            maxo = max(LT.ord(sel));
            if isempty(maxo)||maxo==0, continue; end
            m = nan(1,maxo); se = nan(1,maxo);
            for o = 1:maxo
                v     = LT.VASRating(sel & LT.ord==o);
                m(o)  = mean(v,'omitnan');
                se(o) = std(v,'omitnan')/sqrt(sum(~isnan(v)));
            end
            ls = '-'; if ci==2, ls='--'; end
            errorbar(1:maxo, m, se, [ls 'o'], 'Color', expCols(ei,:), ...
                 'MarkerFaceColor', expCols(ei,:), 'CapSize', 0, ...
                 'HandleVisibility','off');
        end
    end
    xlabel('Conditioning trial within category'); ylabel('Mean VAS \pm SEM');
    box off; title(runs(ri));
    if ri==nRun
        for ei = 1:numel(expList)
            plot(nan,nan,'o-','Color',expCols(ei,:),'MarkerFaceColor',expCols(ei,:), ...
                 'DisplayName',sprintf('Exp %d',expList(ei)));
        end
        plot(nan,nan,'k-','DisplayName','High');
        plot(nan,nan,'k--','DisplayName','Low');
        legend('show','Location','best');
    end
end

% ---------- violin plot means conditioning ----------

cues = {'high','low'};
runs = unique(LT.runid);
nRun = numel(runs);
vw   = 0.35;

figure('Color','w','Name','Conditioning: per-subject means (violins)');
for ri = 1:nRun
    subplot(1,nRun,ri); hold on
    for ci = 1:2
        for ei = 1:numel(expList)
            e    = expList(ei);
            sel  = LT.exp==e & LT.runid==runs(ri) & strcmp(LT.CueAssociation,cues{ci});
            subs = unique(LT.SubID(sel));
            pm   = arrayfun(@(s) mean(LT.VASRating(sel & LT.SubID==s),'omitnan'), subs);
            pm   = pm(~isnan(pm));
            if numel(pm) < 2, continue; end
            xpos = (ci-1)*4 + ei;                       % High:1-3, Low:5-7

            localViolin(xpos, pm, vw, expCols(ei,:));
            jit = (rand(numel(pm),1)-.5)*vw*0.6;
            plot(xpos+jit, pm, '.', 'Color', expCols(ei,:)*.6, ...
                 'MarkerSize', 4, 'HandleVisibility','off');
            plot(xpos, mean(pm), 'ks', 'MarkerFaceColor','k', ...
                 'MarkerSize', 6, 'HandleVisibility','off');
        end
    end
    xticks([2 6]); xticklabels({'High','Low'}); xlim([0 8]);
    ylabel('Per-subject mean VAS'); box off; title(runs(ri));
    if ri==nRun
        for ei = 1:numel(expList)
            fill(nan,nan,expCols(ei,:),'FaceAlpha',.35,'EdgeColor',expCols(ei,:), ...
                 'DisplayName',sprintf('Exp %d',expList(ei)));
        end
        legend('show','Location','best');
    end
end

function localViolin(xc, v, w, col)
    [f,xi] = ksdensity(v);
    f = f / max(f) * w;
    fill([xc+f, fliplr(xc-f)], [xi, fliplr(xi)], col, ...
         'FaceAlpha',.35, 'EdgeColor',col, 'HandleVisibility','off');
end

% ---------- participant variance  ----------
%% per-participant High vs Low with SEM (block rows x experiment cols)
expList = 1:3;
cues    = {'high','low'};
cueCol = [.85 .30 .30; .20 .40 .80];        % High red, Low blue
runs    = unique(LT.runid);                 % ["nocebo","placebo"]
nRun    = numel(runs);

figure('Color','w','Name','Per-participant High vs Low');
p = 0;
for ri = 1:nRun
    for ei = 1:numel(expList)
        e = expList(ei); p = p + 1;
        subplot(nRun, numel(expList), p); hold on

        sel0 = LT.exp==e & LT.runid==runs(ri);
        subs = unique(LT.SubID(sel0));

        % per-subject mean & SEM for each cue
        M = nan(numel(subs),2); S = nan(numel(subs),2);
        for ci = 1:2
            for k = 1:numel(subs)
                v = LT.VASRating(sel0 & LT.SubID==subs(k) & strcmp(LT.CueAssociation,cues{ci}));
                v = v(~isnan(v));
                if ~isempty(v)
                    M(k,ci) = mean(v);
                    S(k,ci) = std(v)/sqrt(numel(v));
                end
            end
        end

        x = 1:numel(subs);                         % subjects in SubID order
        for ci = 1:2
            errorbar(x, M(:,ci), S(:,ci), 'o', 'Color', cueCol(ci,:), ...
                     'MarkerFaceColor', cueCol(ci,:), 'MarkerSize', 3, ...
                     'CapSize', 0, 'LineStyle','none');
        end
        xlim([0 numel(subs)+1]); box off
        title(sprintf('Exp %d — %s', e, runs(ri)));
        if ei==1, ylabel('Mean VAS \pm SEM'); end
        if ri==nRun, xlabel('Subject'); end
        if p==1
            plot(nan,nan,'o','Color',cueCol(1,:),'MarkerFaceColor',cueCol(1,:),'DisplayName','High');
            plot(nan,nan,'o','Color',cueCol(2,:),'MarkerFaceColor',cueCol(2,:),'DisplayName','Low');
            legend('show','Location','best');
        end
    end
end