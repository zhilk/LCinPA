function allResults = run_lmm(behavioralData, studyLabels)
%RUN_LMM  Linear Mixed Model: Placebo vs Nocebo mix trial VAS ratings
%  Runs all three studies sequentially.
%
%  INPUTS:
%    behavioralData  – 3×1 cell array of tables (one per study, already loaded)
%    studyLabels     – 3×1 cell array of study name strings
%
%  OUTPUT:
%    allResults      – 3×1 cell array of result structs (lme, beta, se, t, df, p, d)
%
%  Model per study:
%    VASRating ~ Block + PlaceboFirst + Counterbalancing + TrialC
%               + Block:PlaceboFirst + (1 + Block | SubID)

nStudies = numel(behavioralData);
allResults = cell(nStudies, 1);

for s = 1:nStudies
    allResults{s} = run_study_LMM(behavioralData{s}, studyLabels{s});
end

%% ── Summary table ──────────────────────────────────────────────
fprintf('\n%s\n', repmat('=',1,72));
fprintf('SUMMARY\n');
fprintf('%-25s %7s %7s %8s %8s %8s %8s\n', 'Study','beta','SE','t','df','p','d');
fprintf('%s\n', repmat('-',1,72));
for s = 1:nStudies
    r = allResults{s};
    if r.p<.001, pS='< .001'; else, pS=sprintf('%.3f',r.p); end
    fprintf('%-25s %7.2f %7.2f %8.2f %8.1f %8s %8.3f\n', ...
        studyLabels{s}, r.beta, r.se, r.t, r.df, pS, r.d);
end
fprintf('%s\n', repmat('=',1,72));

end


%% ========================================================================
function results = run_study_LMM(T, studyLabel)

    %% Prepare: test phase, mix trials
    mix = T(strcmp(T.Phase, 'test') & strcmp(T.VisualCategory, 'mix'), :);
    nSub = numel(unique(mix.SubID));

    fprintf('\n%s\n%s — %d subjects, %d mix trials\n%s\n', ...
        repmat('=',1,65), studyLabel, nSub, height(mix), repmat('=',1,65));

    % PlaceboFirst: between-subject recode of BlockOrder
    subs = unique(mix.SubID);
    pfMap = containers.Map('KeyType','double','ValueType','double');
    for i = 1:numel(subs)
        sid = subs(i);
        subRows = mix(mix.SubID == sid & strcmp(mix.Block,'placebo'), :);
        pfMap(sid) = double(subRows.BlockOrder(1) == 1);
    end
    mix.PlaceboFirst = arrayfun(@(x) pfMap(x), mix.SubID);

    fprintf('  PlaceboFirst=1: %d | PlaceboFirst=0: %d\n', ...
        sum(cellfun(@(k) pfMap(k)==1, num2cell(subs))), ...
        sum(cellfun(@(k) pfMap(k)==0, num2cell(subs))));
    fprintf('  Counterbalancing: %s\n', mat2str(unique(mix.Counterbalancing)'));

    % Factors
    mix.SubID_f            = categorical(mix.SubID);
    mix.Block_f            = categorical(mix.Block, {'placebo','nocebo'});
    mix.PlaceboFirstC      = mix.PlaceboFirst - mean(mix.PlaceboFirst);  % mean-centered
    mix.Counterbalancing_f = categorical(mix.Counterbalancing);
    mix.TrialC             = mix.TrialInPhase - mean(mix.TrialInPhase);

    %% Primary model
    formula_full = ['VASRating ~ Block_f + PlaceboFirstC + Counterbalancing_f ' ...
                    '+ TrialC + Block_f:PlaceboFirstC + (1 + Block_f | SubID_f)'];
    formula_red  = ['VASRating ~ Block_f + PlaceboFirstC + Counterbalancing_f ' ...
                    '+ TrialC + Block_f:PlaceboFirstC + (1 | SubID_f)'];
    try
        lme = fitlme(mix, formula_full, 'FitMethod','REML');
        fprintf('\n  Random effects: intercept + slope\n');
    catch
        lme = fitlme(mix, formula_red, 'FitMethod','REML');
        fprintf('\n  Random effects: intercept only (slope singular)\n');
    end

    disp(lme);

    [beta,se,tStat,df,pVal] = extract_block_effect(lme);

    %% ANOVA
    fprintf('── ANOVA ──\n');
    disp(anova(lme));

    %% Cohen's d (subject means)
    sm = grpstats(mix, {'SubID_f','Block_f'}, 'mean', 'DataVars','VASRating');
    pc_m = sm.mean_VASRating(sm.Block_f == 'placebo');
    nc_m = sm.mean_VASRating(sm.Block_f == 'nocebo');
    dv = nc_m - pc_m;
    cohens_d = mean(dv) / std(dv);

    fprintf('── Descriptives (subject means) ──\n');
    fprintf('  Placebo: M = %.2f, SD = %.2f\n', mean(pc_m), std(pc_m));
    fprintf('  Nocebo:  M = %.2f, SD = %.2f\n', mean(nc_m), std(nc_m));
    fprintf('  Cohen''s d = %.3f\n', cohens_d);

    %% Sensitivity 1: Minimal model
    fprintf('\n── Sensitivity: Block only ──\n');
    try
        lme_min = fitlme(mix,'VASRating ~ Block_f + (1+Block_f|SubID_f)','FitMethod','REML');
    catch
        lme_min = fitlme(mix,'VASRating ~ Block_f + (1|SubID_f)','FitMethod','REML');
    end
    [b,s2,t2,d2,p2] = extract_block_effect(lme_min);
    fprintf('  beta=%.2f, SE=%.2f, t(%.1f)=%.2f, p=%.4f\n', b,s2,d2,t2,p2);

    %% Sensitivity 2: Block × Counterbalancing
    fprintf('\n── Sensitivity: Block x Counterbalancing ──\n');
    try
        lme_cb = fitlme(mix, ...
            'VASRating ~ Block_f*Counterbalancing_f + PlaceboFirstC + TrialC + (1+Block_f|SubID_f)', ...
            'FitMethod','REML');
    catch
        lme_cb = fitlme(mix, ...
            'VASRating ~ Block_f*Counterbalancing_f + PlaceboFirstC + TrialC + (1|SubID_f)', ...
            'FitMethod','REML');
    end
    disp(anova(lme_cb));

    %% LRT: full vs null
    fprintf('── LRT: Full vs Null (drop Block) ──\n');
    m_full = fitlme(mix, ...
        ['VASRating ~ Block_f + PlaceboFirstC + Counterbalancing_f + TrialC ' ...
         '+ Block_f:PlaceboFirstC + (1+Block_f|SubID_f)'], 'FitMethod','ML');
    m_null = fitlme(mix, ...
        ['VASRating ~ PlaceboFirstC + Counterbalancing_f + TrialC ' ...
         '+ (1|SubID_f)'], 'FitMethod','ML');
    disp(compare(m_null, m_full));

    %% Diagnostics
    r = residuals(lme);
    fprintf('── Diagnostics ──\n');
    try [~,pN]=lillietest(r); fprintf('  Lilliefors p = %.4f\n',pN);
    catch; try [~,pN]=jbtest(r); fprintf('  Jarque-Bera p = %.4f\n',pN);
    catch; end; end
    fprintf('  Residual SD = %.2f\n', std(r));

    %% Store
    results.lme=lme; results.beta=beta; results.se=se;
    results.t=tStat; results.df=df; results.p=pVal; results.d=cohens_d;
end


%% ========================================================================
%  HELPERS
%% ========================================================================
function [beta,se,tStat,df,pVal] = extract_block_effect(lme)
    [~,~,stats] = fixedEffects(lme, 'DFMethod','Satterthwaite');
    idx = find(contains(stats.Name,'nocebo') & ~contains(stats.Name,':'), 1);
    beta=stats.Estimate(idx); se=stats.SE(idx); tStat=stats.tStat(idx);
    df=stats.DF(idx); pVal=stats.pValue(idx);
end
