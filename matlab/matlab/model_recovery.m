function model_recovery(nRep)
%MODEL_RECOVERY  Model recovery (confusion matrices) for the six-model set.
%
%  Per generating model and subject: nRep (default 3) synthetic datasets
%      CR_sim = b0 + b1*R(theta_ML) + N(0, sigma_ML^2)   (all trials)
%  at the subject's own ML parameters on the subject's own trial sequence.
%  All six models are refit to each dataset (valid-test scoring); two
%  winners are recorded:
%    best_bic    lowest BIC on valid test trials;
%    best_morph  highest summed log likelihood on valid morph trials,
%                evaluated at each model's valid-test ML parameters
%                (parameters are not refit to the morph trials).
%  Ties resolve in model order null, rw, bcc, rwbcc, simple, percept.
%
%  Output: results/recovery_model_confusion.csv (one row per dataset).
%  Seed rng(11); plain loops (subject > generator > replicate), exactly
%  reproducible.

    if nargin < 1, nRep = 3; end
    root = fileparts(fileparts(mfilename('fullpath')));
    subs = load_subjects(root);
    g = model_grids();
    models = {'null', 'rw', 'bcc', 'rwbcc', 'simple', 'percept'};
    nM = numel(models);

    rng(11, 'twister');
    rows = [];
    for i = 1:numel(subs)
        sub = subs{i};
        vm = sub.valid & sub.isMix;
        P = cell(1, nM); PP = cell(1, nM); F = cell(1, nM);
        for mi = 1:nM
            [P{mi}, PP{mi}] = predict_model(sub, models{mi}, g);
            F{mi} = fit_on_grid(P{mi}, PP{mi}, sub.CR, sub.vt, g.npar.(models{mi}));
        end
        for gi = 1:nM
            f = F{gi};
            base = f.b0 + f.b1*P{gi}(:, f.idx);
            for rep = 1:nRep
                CRsim = base + f.sigma*randn(sub.n, 1);
                bic = zeros(1, nM); mll = zeros(1, nM);
                for mi = 1:nM
                    r = fit_on_grid(P{mi}, PP{mi}, CRsim, sub.vt, g.npar.(models{mi}));
                    bic(mi) = r.BIC;
                    pr = r.b0 + r.b1*P{mi}(:, r.idx);
                    rm = CRsim(vm) - pr(vm);
                    mll(mi) = sum(-0.5*log(2*pi) - log(r.sigma) ...
                                  - 0.5*(rm/r.sigma).^2);
                end
                [~, wb] = min(bic);
                [~, wm] = max(mll);
                row = struct('generator', string(models{gi}), ...
                             'study', sub.study, 'SubID', sub.SubID, ...
                             'rep', rep, 'best_bic', string(models{wb}), ...
                             'best_morph', string(models{wm}));
                if isempty(rows), rows = row; else, rows(end+1) = row; end %#ok<AGROW>
            end
        end
        if mod(i, 10) == 0
            fprintf('model_recovery %d/%d\n', i, numel(subs));
        end
    end
    out = fullfile(root, 'results', 'recovery_model_confusion.csv');
    writetable(struct2table(rows), out);
    fprintf('written %s\n', out);
end
