%% ========================================================================
%  Model diagnostics: kappa identifiability, observation-model parameters,
%  an example subject, and rating bias by stimulus type and block.
%% ========================================================================

s = 1;                                   % study index for panels 1–3
R = allStudyResults{s};
ok = ~isnan(R.BIC_RW);
idx = find(ok);
n   = numel(idx);


% %% ── Example subject: ratings and RW weights by block ──────────────
% jx  = randi(numel(idx));
% ix  = idx(jx);
% sid = R.subs(ix);
% rRW = R.allRW{ix};
% 
% sTx = behavioralData{s}(behavioralData{s}.SubID==sid,:);
% sTx = sTx(sTx.VASResponse==1 & ~isnan(sTx.VASRating),:);
% sTx = sortrows(sTx,'TrialGlobal');
% 
% W = rRW.w_t;
% fprintf('sub %d: W has %d rows, table has %d trials\n', sid, size(W,1), height(sTx));
% 
% if size(W,1) < height(sTx)
%     sTx = sTx(~strcmp(sTx.Phase,'conditioning'), :);   % model ran test-only
%     fprintf('  -> subset to %d test trials\n', height(sTx));
% end
% 
% blk = string(sTx.Block);
% ub  = unique(blk,'stable');
% isC = strcmp(sTx.Phase,'conditioning');
% US  = sTx.TargetVAS/100;
% CR  = sTx.VASRating/100;
% isM = abs(sTx.x_face-0.5)<1e-9;
% 
% figure('Color','w','Position',[100 100 1250 420]);
% for b = 1:numel(ub)
%     m  = blk==ub(b);
%     tb = (1:nnz(m))';
%     Ub = US(m); Cb = CR(m); Mb = isM(m); Wb = W(m,:); Cc = isC(m);
% 
%     lvl  = unique(Ub(~Mb));
%     loU  = min(lvl); hiU = max(lvl);
%     isLo = ~Mb & abs(Ub-loU)<1e-9;
%     isHi = ~Mb & abs(Ub-hiU)<1e-9;
% 
%     subplot(1,numel(ub),b); hold on
%     yline(loU,'Color',[0 .6 .2],'LineStyle',':');
%     yline(hiU,'Color',[.85 .1 .1],'LineStyle',':');
%     yline(0.5,'Color',[.1 .3 .9],'LineStyle',':');
%     plot(tb, Wb(:,1), '-', 'Color',[.3 .3 .3], 'LineWidth',1.4);
%     plot(tb, Wb(:,2), '-', 'Color',[.65 .65 .65],'LineWidth',1.4);
%     plot(tb(isLo), Cb(isLo), '.', 'Color',[0 .6 .2],  'MarkerSize',13);
%     plot(tb(isHi), Cb(isHi), '.', 'Color',[.85 .1 .1],'MarkerSize',13);
%     plot(tb(Mb),   Cb(Mb),   '.', 'Color',[.1 .3 .9], 'MarkerSize',13);
%     if any(Cc); xline(nnz(Cc)+0.5,'k--'); end
%     ylim([0 1]); box off
%     xlabel('trial in block'); ylabel('normalized (0–1)'); title(ub(b))
%     if b==1
%         legend({'low US','high US','morph US','w_{face}','w_{house}', ...
%                 'CR low','CR high','CR morph'},'Location','best','FontSize',7);
%     end
% end
% sgtitle(sprintf('sub %d  |  winner: %s  |  \\eta=%.3f  \\kappa_{BCC}=%.2f  \\alpha=%.2f', ...
%     sid, R.bestModel{ix}, R.eta_RW(ix), R.kappa_BCC(ix), R.alpha_LCM(ix)))
% 
% 
% %% ── Rating bias by stimulus type and block, all studies ───────────
% nS   = numel(allStudyResults);
% Dall = cell(nS,1);
% cols = [0 .6 .2; .85 .1 .1; .1 .3 .9];      % low, high, morph
% 
% figure('Color','w','Position',[80 80 1400 420]);
% for ss = 1:nS
%     Rs  = allStudyResults{ss};
%     okS = ~isnan(Rs.BIC_RW); idxS = find(okS);
%     D   = nan(numel(idxS), 6);              % [lo hi morph] x [placebo nocebo]
% 
%     for j = 1:numel(idxS)
%         sid = Rs.subs(idxS(j));
%         sT  = behavioralData{ss}(behavioralData{ss}.SubID==sid,:);
%         sT  = sT(sT.VASResponse==1 & ~isnan(sT.VASRating),:);
%         US  = sT.TargetVAS; CR = sT.VASRating;
%         blk = string(sT.Block);
%         isM = abs(sT.x_face-0.5)<1e-9;
%         c = 0;
%         for bb = ["placebo","nocebo"]
%             m   = blk==bb;
%             lvl = unique(US(m & ~isM));
%             for typ = 1:3
%                 c = c+1;
%                 switch typ
%                     case 1, sel = m & ~isM & US==min(lvl);
%                     case 2, sel = m & ~isM & US==max(lvl);
%                     case 3, sel = m & isM;
%                 end
%                 if any(sel); D(j,c) = mean(CR(sel) - US(sel)); end
%             end
%         end
%     end
%     Dall{ss} = D;
% 
%     subplot(1,nS,ss); hold on
%     for c = 1:6
%         x = D(~isnan(D(:,c)), c);
%         if numel(x) < 2; continue; end
%         gy = linspace(min(x), max(x), 100);
%         try, f = ksdensity(x, gy); catch, [f,gy] = hist(x,12); end %#ok
%         f   = 0.35 * f / max(f);
%         col = cols(mod(c-1,3)+1, :);
%         fill([c+f, fliplr(c-f)], [gy, fliplr(gy)], col, ...
%              'FaceAlpha',0.35, 'EdgeColor',col);
%         plot([c-.35 c+.35], [median(x) median(x)], 'k-','LineWidth',1.5);
%         plot(c + 0.06*randn(numel(x),1), x, 'k.','MarkerSize',4);
%     end
%     yline(0,'k--'); xline(3.5,'k:');
%     set(gca,'XTick',1:6,'XTickLabel',{'lo','hi','mix','lo','hi','mix'});
%     box off; ylabel('CR - US (VAS)');
%     title(sprintf('%s (n=%d)', studyLabels{ss}, numel(idxS)))
%     text(2, max(ylim)*0.95, 'placebo', 'HorizontalAlignment','center');
%     text(5, max(ylim)*0.95, 'nocebo',  'HorizontalAlignment','center');
% end
% sgtitle('Rating bias relative to delivered pain, by stimulus and block')


%% ── kappa across models, per participant ─────────────────────────
s  = 1;
R  = allStudyResults{s};
ok = ~isnan(R.BIC_BCC);
idx = find(ok);
n   = numel(idx);
tol = 1e-4;

% candidate models and where their per-subject kappa lives
cand = { 'BCC',    'allBCC'
         'RW+BCC', 'allHYB'
         'RW',     'allRW'
         'LCM',    'allLCM' };

K = []; L = {};
for c = 1:size(cand,1)
    fld = cand{c,2};
    if ~isfield(R, fld); continue; end
    v = nan(n,1); has = false;
    for j = 1:n
        r = R.(fld){idx(j)};
        if ~isempty(r) && isfield(r,'kappa'); v(j) = r.kappa; has = true; end
    end
    if has; K(:,end+1) = v; L{end+1} = cand{c,1}; end %#ok
end

figure('Color','w','Position',[100 100 1000 430]);

% (1) per-participant kappa, paired
subplot(1,2,1); hold on
for i = 1:n
    plot(1:numel(L), K(i,:), '-', 'Color',[.85 .85 .85]);
end
boxplot(K, 'Labels',L, 'Colors','k', 'Symbol','', 'Widths',0.5);
plot(1:numel(L), median(K,1,'omitnan'), 'rd','MarkerFaceColor','r','MarkerSize',8);
ylim([-0.05 1.05]); box off
ylabel('\kappa (weight on expectation)')
title(sprintf('%s — \\kappa per participant (n=%d)', studyLabels{s}, n))

% (2) boundary solutions
subplot(1,2,2);
bnd = [mean(K<=tol,1); mean(K>=1-tol,1)]' * 100;
bar(bnd,'stacked'); set(gca,'XTickLabel',L);
ylabel('% subjects at bound'); box off
legend({'\kappa\approx0','\kappa\approx1'},'Location','best');
title('boundary solutions')

for c = 1:numel(L)
    fprintf('%-8s median %.2f   %.0f%% at 0   %.0f%% at 1\n', ...
        L{c}, median(K(:,c),'omitnan'), ...
        100*mean(K(:,c)<=tol), 100*mean(K(:,c)>=1-tol));
end

%% ── Cause vectors, one subject ───────────────────────────────────
s   = 1;
R   = allStudyResults{s};
ok  = ~isnan(R.BIC_LCM); idx = find(ok);
ix  = idx(randi(numel(idx)));          % or pick: ix = find(R.subs==8);
sid = R.subs(ix);
r   = R.allLCM{ix};

mu = (r.SumF + r.a*r.mu0) ./ (r.Nk(:) + r.a);   % [K x 3]
K  = size(mu,1);

C = table((1:K)', r.Nk(:), mu(:,1), mu(:,2), mu(:,3), ...
          'VariableNames', {'k','Nk','muUS','muFace','muHouse'});
C.fromCond = ismember(C.k, r.condCauses);
C.type = strings(K,1);
C.type(mu(:,2)>0.65) = "face";
C.type(mu(:,3)>0.65) = "house";
C.type(abs(mu(:,2)-0.5)<0.15 & abs(mu(:,3)-0.5)<0.15) = "morph";
C.type(C.type=="") = "mixed";

fprintf('\nsubject %d — %d causes\n', sid, K);
disp(C);

% cause vectors in CS x US space
figure('Color','w','Position',[100 100 560 480]); hold on
cols = containers.Map({'face','house','morph','mixed'}, ...
                      {[0 .3 .9],[0 .6 .2],[.85 .1 .1],[.5 .5 .5]});
for k = 1:K
    sz = 20 + 6*C.Nk(k);
    scatter(mu(k,2), mu(k,1), sz, cols(char(C.type(k))), 'filled', ...
            'MarkerEdgeColor','k', 'LineWidth', 0.5+1.0*C.fromCond(k));
    text(mu(k,2)+0.02, mu(k,1), sprintf('%d (n=%d)', k, C.Nk(k)), 'FontSize',8);
end
xline(0.5,'k:'); yline(0.5,'k:');
xlim([-0.1 1.1]); ylim([0 1]); box off
xlabel('\mu_{face}  (0 = house, 1 = face, 0.5 = morph)')
ylabel('\mu_{US}  (learned pain)')
title(sprintf('sub %d — cause vectors (size \\propto N_k, thick edge = from conditioning)', sid))

%% ── Number of causes per subject ─────────────────────────────────
s  = 1;
R  = allStudyResults{s};
ok = ~isnan(R.BIC_LCM);
idx = find(ok);
n   = numel(idx);

nCause = nan(n,1); nMorph = nan(n,1); nCS = nan(n,1);
for j = 1:n
    r = R.allLCM{idx(j)};
    if isempty(r); continue; end
    mu = (r.SumF + r.a*r.mu0) ./ (r.Nk(:) + r.a);
    nCause(j) = size(mu,1);
    isMorph   = abs(mu(:,2)-0.5)<0.15 & abs(mu(:,3)-0.5)<0.15;
    nMorph(j) = nnz(isMorph);
    nCS(j)    = nnz(~isMorph);
end

figure('Color','w','Position',[100 100 1150 400]);

% (1) total causes
subplot(1,3,1);
histogram(nCause, 0.5:1:max(nCause)+0.5, 'FaceColor',[.4 .5 .8]); box off
xline(6,'r--','expected (3 x 2 blocks)');
xlabel('# causes'); ylabel('# subjects')
title(sprintf('total: median %d, range %d–%d', ...
    median(nCause,'omitnan'), min(nCause), max(nCause)))

% (2) split CS vs morph causes
subplot(1,3,2);
bar([histcounts(nCS, 0.5:1:10.5); histcounts(nMorph, 0.5:1:10.5)]');
set(gca,'XTick',1:10); box off
xlabel('# causes of that kind'); ylabel('# subjects')
legend({'face/house causes','morph causes'},'Location','best')
title('cause types')

% (3) per subject
subplot(1,3,3);
bar(R.subs(ok), [nCS, nMorph], 'stacked'); box off
yline(6,'r--');
xlabel('subject'); ylabel('# causes')
legend({'face/house','morph'},'Location','best')
title('causes per subject')

sgtitle(sprintf('%s — latent causes opened (n=%d, structure at fixed \\alpha=1)', ...
    studyLabels{s}, n))

fprintf('causes: median %.1f (range %d–%d)\n', median(nCause,'omitnan'), ...
    min(nCause), max(nCause));
fprintf('%.0f%% of subjects form exactly 6 (the canonical 3 x 2)\n', 100*mean(nCause==6));
fprintf('%.0f%% at K_cap (10)\n', 100*mean(nCause>=10));