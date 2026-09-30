s = 1;
R = allStudyResults{s};
ok = ~isnan(R.BIC_RW);
idx = find(ok);
n   = numel(idx);

kGrid = linspace(0,1,51);
LLk = nan(n, numel(kGrid), 3);
b0  = nan(n,3); b1 = nan(n,3); sg = nan(n,3);

for j = 1:n
    i   = idx(j);
    sid = R.subs(i);
    sT  = behavioralData{s}(behavioralData{s}.SubID==sid,:);
    sT  = sT(sT.VASResponse==1 & ~isnan(sT.VASRating),:);
    sT  = sortrows(sT,'TrialGlobal');
    US  = sT.TargetVAS/100;  CR = sT.VASRating;

    mods = {R.allRW{i}, R.allBCC{i}, R.allLCM{i}};
    for m = 1:3
        r = mods{m};
        if isempty(r); continue; end
        b0(j,m) = r.b0; b1(j,m) = r.b1; sg(j,m) = r.sigma;
        if ~isfield(r,'V'); continue; end
        V = r.V(:);
        if numel(V)~=numel(US); continue; end
        for q = 1:numel(kGrid)
            LLk(j,q,m) = rescaled_LL(kGrid(q)*V + (1-kGrid(q))*US, CR);
        end
    end
end

L = {'RW','BCC','LCM'};
figure('Color','w','Position',[100 100 1250 800]);

% ── row 1: kappa likelihood profiles ──
for m = 1:3
    subplot(2,3,m); hold on
    Lm   = LLk(:,:,m);
    Lrel = Lm - max(Lm,[],2);
    plot(kGrid, Lrel', '-', 'Color',[.75 .75 .75 .5]);
    plot(kGrid, median(Lrel,1,'omitnan'), 'k-', 'LineWidth',2);
    yline(-1.92,'r--');
    ylim([-15 1]); box off
    xlabel('\kappa'); ylabel('\DeltaLL from subject max')
    title(sprintf('%s — median drop at \\kappa=0: %.1f', ...
        L{m}, median(Lrel(:,1),'omitnan')))
end

%%
% ── row 2: observation-model parameters (outliers hidden) ──
subplot(2,3,4);
boxplot(b0,'Labels',L,'Symbol',''); yline(0,'k:'); box off
ylim([min(prctile(b0,2)), max(prctile(b0,98))])
ylabel('b_0 (intercept, VAS)'); title('rescaling offset (2–98 pct)')

subplot(2,3,5);
boxplot(b1,'Labels',L,'Symbol',''); yline(0,'r:'); box off
ylim([min(prctile(b1,2)), max(prctile(b1,98))])
ylabel('b_1 (slope)'); title('rescaling gain (2–98 pct)')

subplot(2,3,6);
boxplot(sg,'Labels',L,'Symbol',''); box off
ylim([min(prctile(sg,2)), max(prctile(sg,98))])
ylabel('\sigma (residual SD, VAS)'); title('unexplained noise (2–98 pct)')

sgtitle(sprintf('%s — \\kappa identifiability and observation model (n=%d)', ...
    studyLabels{s}, n))

fprintf('median b0: %s\n', sprintf('%7.1f', median(b0,1,'omitnan')));
fprintf('median b1: %s\n', sprintf('%7.1f', median(b1,1,'omitnan')));
fprintf('median sg: %s\n', sprintf('%7.1f', median(sg,1,'omitnan')));
fprintf('extreme b1 (|b1|>500): %d subjects\n', nnz(any(abs(b1)>500,2)));

s = 3;
R = allStudyResults{s};
outDir = fullfile(pwd, 'plots');                  % existing subfolder

ok  = ~isnan(R.BIC_RW);
idx = find(ok);

for j = 1:numel(idx)
    i   = idx(j);
    sid = R.subs(i);
    rRW = R.allRW{i};

    sT = behavioralData{s}(behavioralData{s}.SubID==sid,:);
    sT = sT(sT.VASResponse==1 & ~isnan(sT.VASRating),:);
    sT = sortrows(sT,'TrialGlobal');

    W   = rRW.w_t;
    blk = string(sT.Block);
    ub  = unique(blk,'stable');
    isC = strcmp(sT.Phase,'conditioning');
    US  = sT.TargetVAS/100;
    CR  = sT.VASRating/100;
    isM = abs(sT.x_face-0.5)<1e-9;

    f = figure('Color','w','Position',[100 100 1250 420],'Visible','off');

    for b = 1:numel(ub)
        m  = blk==ub(b);
        tb = (1:nnz(m))';
        Ub = US(m); Cb = CR(m); Mb = isM(m); Wb = W(m,:); Cc = isC(m);

        lvl  = unique(Ub(~Mb));
        loU  = min(lvl); hiU = max(lvl);
        isLo = ~Mb & abs(Ub-loU)<1e-9;
        isHi = ~Mb & abs(Ub-hiU)<1e-9;

        subplot(1,numel(ub),b); hold on
        yline(loU,'Color',[0 .6 .2],'LineStyle',':');
        yline(hiU,'Color',[.85 .1 .1],'LineStyle',':');
        yline(0.5,'Color',[.1 .3 .9],'LineStyle',':');
        plot(tb, Wb(:,1), '-', 'Color',[.3 .3 .3], 'LineWidth',1.4);
        plot(tb, Wb(:,2), '-', 'Color',[.65 .65 .65],'LineWidth',1.4);
        plot(tb(isLo), Cb(isLo), '.', 'Color',[0 .6 .2],  'MarkerSize',13);
        plot(tb(isHi), Cb(isHi), '.', 'Color',[.85 .1 .1],'MarkerSize',13);
        plot(tb(Mb),   Cb(Mb),   '.', 'Color',[.1 .3 .9], 'MarkerSize',13);
        xline(nnz(Cc)+0.5,'k--');
        ylim([0 1]); box off
        xlabel('trial in block'); ylabel('normalized (0–1)'); title(ub(b))
        if b==1
            legend({'low US','high US','morph US','w_{face}','w_{house}', ...
                    'CR low','CR high','CR morph'},'Location','best','FontSize',7);
        end
    end

    sgtitle(sprintf(['sub %d  |  winner: %s  |  \\eta=%.3f  ' ...
                     '\\kappa_{RW}=%.2f  \\kappa_{BCC}=%.2f  \\kappa_{LCM}=%.2f  \\alpha=%.2f'], ...
        sid, R.bestModel{i}, R.eta_RW(i), R.kappa_RW(i), ...
        R.kappa_BCC(i), R.kappa_LCM(i), R.alpha_LCM(i)))

    exportgraphics(f, fullfile(outDir, sprintf('study%d_sub%03d.png', s, sid)), ...
                   'Resolution',150);
    close(f);
end
fprintf('saved %d figures to %s\n', numel(idx), outDir);