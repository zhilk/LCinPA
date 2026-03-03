function results = lcm_infer(X, opts)
%LCM_INFER  Local MAP inference for latent cause model.
%
%  Based on Gershman et al. (2010, 2017) and the sjgershm/LCM repository.
%  Adapted for the LCPA placebo/nocebo paradigm.
%
%  INPUTS:
%    X    – [T x D] stimulus matrix. Column 1 = US (normalized 0-1),
%           columns 2:D = CS features (face, house).
%    opts – struct with fields:
%           .alpha      – CRP concentration parameter (default 1)
%           .stickiness – stickiness for most recent cause (default 0)
%           .lambda     – temporal decay rate for cause counts (default 0)
%                         Nk = exp(-lambda) * Nk on each trial (Heesink et al., 2024)
%           .a          – Beta prior pseudocount (default 1)
%           .b          – Beta prior pseudocount (default 1)
%           .K          – initial max number of causes (default 10)
%           .sigma_cs   – CS feature variance (default 0.5)
%           .sigma_us   – US feature variance (default 0.5)
%
%  OUTPUTS:
%    results.V    – [T x 1] US predictions (from CS-only posterior)
%    results.post – [T x K] posterior cause probabilities
%    results.z    – [T x 1] MAP cause assignments
%    results.K    – final number of active causes
%    results.opts – options used

    if nargin < 2; opts = struct(); end
    opts = set_defaults(opts);
    results.opts = opts;

    [T, D] = size(X);
    K = opts.K;
    if opts.alpha == 0; K = 1; end

    % Sufficient statistics per cause
    Nk    = zeros(1, K);          % observation count per cause
    SumF  = zeros(K, D);          % running sum of features per cause
    SumF2 = zeros(K, D);          % running sum of squared features

    % Feature-specific variance
    sigma2 = ones(1, D);
    sigma2(1)   = opts.sigma_us;  % US variance
    sigma2(2:D) = opts.sigma_cs;  % CS variance

    % Prior mean for all features
    mu0 = 0.5;
    a_pseudo = opts.a;            % pseudocount for prior

    % Storage
    results.V    = zeros(T, 1);
    results.post = zeros(T, K);
    results.z    = zeros(T, 1);
    z_prev = 1;                   % most recently active cause

    for t = 1:T
        xt = X(t, :);             % [1 x D]: US, CS1, CS2

        % --- Temporal decay of cause counts (Heesink et al., 2024) ---
        if opts.lambda > 0
            Nk = exp(-opts.lambda) * Nk;
        end

        if opts.alpha > 0
            % --- CRP prior ---
            prior = Nk;
            prior(z_prev) = prior(z_prev) + opts.stickiness;
            idx_new = find(Nk == 0, 1);
            if isempty(idx_new)
                % Expand capacity
                K = K + 1;
                Nk(K) = 0;
                SumF(K,:) = 0;
                SumF2(K,:) = 0;
                results.post(:,K) = 0;
                prior(K) = 0;
                idx_new = K;
            end
            prior(idx_new) = opts.alpha;
            prior = prior / sum(prior);

            % --- Likelihood per cause per feature (Gaussian) ---
            lik = zeros(1, K);
            lik_cs = zeros(1, K);
            for k = 1:K
                log_lik_full = 0;
                log_lik_cs   = 0;
                for d = 1:D
                    if Nk(k) > 0
                        mu_hat = (SumF(k,d) + a_pseudo * mu0) / (Nk(k) + a_pseudo);
                    else
                        mu_hat = mu0;
                    end
                    var_pred = sigma2(d) + sigma2(d) / (Nk(k) + a_pseudo);
                    ll_d = -0.5 * log(2*pi*var_pred) - 0.5 * (xt(d) - mu_hat)^2 / var_pred;
                    log_lik_full = log_lik_full + ll_d;
                    if d >= 2
                        log_lik_cs = log_lik_cs + ll_d;
                    end
                end
                lik(k)    = exp(log_lik_full);
                lik_cs(k) = exp(log_lik_cs);
            end

            % --- CS-only posterior (for prediction) ---
            post_cs = prior .* lik_cs;
            post_cs = post_cs / max(sum(post_cs), eps);

            % --- Full posterior (CS + US, for cause assignment) ---
            post_full = prior .* lik;
            post_full = post_full / max(sum(post_full), eps);
        else
            % alpha = 0: single cause, no inference needed
            post_cs   = zeros(1, K);  post_cs(1)   = 1;
            post_full = zeros(1, K);  post_full(1)  = 1;
        end

        % --- US prediction: posterior-weighted average of cause means ---
        V_t = 0;
        for k = 1:K
            if Nk(k) > 0
                mu_us_k = (SumF(k,1) + a_pseudo * mu0) / (Nk(k) + a_pseudo);
            else
                mu_us_k = mu0;
            end
            V_t = V_t + post_cs(k) * mu_us_k;
        end
        results.V(t) = V_t;

        results.post(t, 1:K) = post_full;

        % --- MAP assignment ---
        [~, z_map] = max(post_full);
        results.z(t) = z_map;
        z_prev = z_map;

        % --- Update sufficient statistics ---
        Nk(z_map)       = Nk(z_map) + 1;
        SumF(z_map, :)  = SumF(z_map, :) + xt;
        SumF2(z_map, :) = SumF2(z_map, :) + xt.^2;
    end

    % Trim unused causes
    active = sum(results.post, 1) > 0;
    results.post = results.post(:, active);
    results.K = sum(active);

    % Store final sufficient statistics (for cause assignment analysis)
    results.Nk   = Nk(active);
    results.SumF = SumF(active, :);
    results.mu0  = mu0;
    results.a_pseudo = a_pseudo;
end


%% ========================================================================
function opts = set_defaults(opts)
    if ~isfield(opts,'alpha')      || isempty(opts.alpha);      opts.alpha = 1;     end
    if ~isfield(opts,'stickiness') || isempty(opts.stickiness); opts.stickiness = 0; end
    if ~isfield(opts,'lambda')     || isempty(opts.lambda);     opts.lambda = 0;    end
    if ~isfield(opts,'a')          || isempty(opts.a);          opts.a = 1;         end
    if ~isfield(opts,'b')          || isempty(opts.b);          opts.b = 1;         end
    if ~isfield(opts,'K')          || isempty(opts.K);          opts.K = 10;        end
    if ~isfield(opts,'sigma_cs')   || isempty(opts.sigma_cs);   opts.sigma_cs = 0.5; end
    if ~isfield(opts,'sigma_us')   || isempty(opts.sigma_us);   opts.sigma_us = 0.5; end
end
