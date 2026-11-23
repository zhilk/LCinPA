function fit_all_subjects()
%FIT_ALL_SUBJECTS  Fit all six models to every subject (N = 156).
%
%  Valid-test scoring, grids of record. Output results/matlab_fits.csv,
%  one row per subject x model:
%    par1 / par2:  rw eta / -;  bcc kappa / -;  rwbcc eta / kappa;
%                  simple sigma_us / -;  percept alpha / sigma_us
%    effPred: model-predicted morph effect (VAS, nocebo minus placebo)
%    effObs:  observed effect on valid morph trials (same for all models)

    root = fileparts(fileparts(mfilename('fullpath')));
    subs = load_subjects(root);
    models = {'null', 'rw', 'bcc', 'rwbcc', 'simple', 'percept'};

    rows = [];
    for i = 1:numel(subs)
        sub = subs{i};
        effObs = mean(sub.CR(sub.isMix & sub.valid & sub.blocks == "nocebo")) ...
               - mean(sub.CR(sub.isMix & sub.valid & sub.blocks == "placebo"));
        for mi = 1:numel(models)
            m = models{mi};
            switch m
                case 'null',    res = fit_null(sub.T);
                case 'rw',      res = fit_rw(sub.T);
                case 'bcc',     res = fit_bcc(sub.T);
                case 'rwbcc',   res = fit_rwbcc(sub.T);
                case 'simple',  res = fit_lcm_simple(sub.T);
                case 'percept', res = fit_lcm_percept(sub.T);
            end
            switch m
                case 'null',    p = [NaN, NaN];
                case 'rw',      p = [res.eta, NaN];
                case 'bcc',     p = [res.kappa, NaN];
                case 'rwbcc',   p = [res.eta, res.kappa];
                case 'simple',  p = [res.sigma_us, NaN];
                case 'percept', p = [res.alpha, res.sigma_us];
            end
            r = struct('study', sub.study, 'SubID', sub.SubID, ...
                       'model', string(m), 'par1', p(1), 'par2', p(2), ...
                       'LL', res.LL, 'BIC', res.BIC, ...
                       'b0', res.b0, 'b1', res.b1, 'sigma', res.sigma, ...
                       'effPred', eff_pred(sub, res.R, res.b0, res.b1), ...
                       'effObs', effObs);
            if isempty(rows), rows = r; else, rows(end+1) = r; end %#ok<AGROW>
        end
        if mod(i, 20) == 0
            fprintf('fit_all_subjects %d/%d\n', i, numel(subs));
        end
    end
    out = fullfile(root, 'results', 'matlab_fits.csv');
    writetable(struct2table(rows), out);
    fprintf('written %s\n', out);
end
