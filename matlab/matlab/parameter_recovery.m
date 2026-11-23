function parameter_recovery(nRep)
%PARAMETER_RECOVERY  Parameter recovery for all models with free model
%  parameters (rw, bcc, rwbcc, simple, percept).
%
%  Per subject and model: nRep (default 10) synthetic datasets are
%  generated from the subject's own ML parameters on the subject's own
%  trial sequence,
%      CR_sim = b0 + b1*R(theta_ML) + N(0, sigma_ML^2)   (all trials),
%  and refit with the identical grid pipeline (valid-test scoring).
%  Because predictions depend on the stimuli only, each grid point's
%  forward pass is computed once per subject (predict_model) and every
%  refit is the closed-form rating fit (fit_on_grid).
%
%  Recorded per synthetic dataset: true vs recovered parameters; for the
%  two latent-cause models additionally the true vs recovered predicted
%  morph effect, and for the simple model the assignment phenotype.
%
%  Outputs (results/): recovery_param_simple.csv, recovery_param_percept.csv,
%  recovery_param_rwbcc.csv (models rw, bcc, rwbcc in one file).
%  Seed rng(7); plain loops (subject > model > replicate), so the run is
%  exactly reproducible.

    if nargin < 1, nRep = 10; end
    root = fileparts(fileparts(mfilename('fullpath')));
    subs = load_subjects(root);
    g = model_grids();
    models = {'rw', 'bcc', 'rwbcc', 'simple', 'percept'};

    rng(7, 'twister');
    rs = []; rp = []; rb = [];
    for i = 1:numel(subs)
        sub = subs{i};
        for mi = 1:numel(models)
            m = models{mi};
            [pred, params] = predict_model(sub, m, g);
            f = fit_on_grid(pred, params, sub.CR, sub.vt, g.npar.(m));
            base = f.b0 + f.b1*pred(:, f.idx);
            for rep = 1:nRep
                CRsim = base + f.sigma*randn(sub.n, 1);
                r = fit_on_grid(pred, params, CRsim, sub.vt, g.npar.(m));
                switch m
                    case 'simple'
                        row = struct('study', sub.study, 'SubID', sub.SubID, ...
                            'rep', rep, 'true_sus', f.p(1), 'rec_sus', r.p(1), ...
                            'true_eff', eff_pred(sub, pred(:, f.idx), f.b0, f.b1), ...
                            'rec_eff',  eff_pred(sub, pred(:, r.idx), r.b0, r.b1), ...
                            'true_ph', phenotype_simple(sub, f.p(1)), ...
                            'rec_ph',  phenotype_simple(sub, r.p(1)));
                        if isempty(rs), rs = row; else, rs(end+1) = row; end %#ok<AGROW>
                    case 'percept'
                        row = struct('study', sub.study, 'SubID', sub.SubID, ...
                            'rep', rep, ...
                            'true_alpha', f.p(1), 'rec_alpha', r.p(1), ...
                            'true_sus',   f.p(2), 'rec_sus',   r.p(2), ...
                            'true_eff', eff_pred(sub, pred(:, f.idx), f.b0, f.b1), ...
                            'rec_eff',  eff_pred(sub, pred(:, r.idx), r.b0, r.b1));
                        if isempty(rp), rp = row; else, rp(end+1) = row; end %#ok<AGROW>
                    otherwise   % rw / bcc / rwbcc
                        switch m
                            case 'rw',    tp = [f.p(1) NaN]; rr = [r.p(1) NaN];
                            case 'bcc',   tp = [NaN f.p(1)]; rr = [NaN r.p(1)];
                            case 'rwbcc', tp = f.p;          rr = r.p;
                        end
                        row = struct('study', sub.study, 'SubID', sub.SubID, ...
                            'model', string(m), 'rep', rep, ...
                            'true_eta',   tp(1), 'rec_eta',   rr(1), ...
                            'true_kappa', tp(2), 'rec_kappa', rr(2));
                        if isempty(rb), rb = row; else, rb(end+1) = row; end %#ok<AGROW>
                end
            end
        end
        if mod(i, 10) == 0
            fprintf('parameter_recovery %d/%d\n', i, numel(subs));
        end
    end
    writetable(struct2table(rs), fullfile(root, 'results', 'recovery_param_simple.csv'));
    writetable(struct2table(rp), fullfile(root, 'results', 'recovery_param_percept.csv'));
    writetable(struct2table(rb), fullfile(root, 'results', 'recovery_param_rwbcc.csv'));
    fprintf('parameter recovery done (%d subjects x %d reps)\n', numel(subs), nRep);
end
