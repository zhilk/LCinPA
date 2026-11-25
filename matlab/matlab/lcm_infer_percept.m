function results = lcm_infer_percept(X, opts)

%  Local MAP inference for the CRP latent-cause model, one block of trials.
%  opts.sigma_cs / .sigma_us are standard deviations (squared internally).
%
%  Each trial t with features x_t = [US, x_face]:
%    1. prior over causes: occupied causes proportional to their visit
%       counts Nk (plus opts.stickiness on the previous cause); opts.alpha
%       on the first empty slot (capacity opts.K; the last slot is reused
%       when capacity is full). alpha = 0 reduces to a single cause.
%    2. Gaussian feature likelihood per cause: means
%       (SumF + a*mu0)/(Nk + a), variances sigma^2 * (1 + 1/(Nk + a)).
%    3. read-outs: V(t) = sum_k post_cs(k) * mu_US(k) from the CS-only
%       posterior (expectation before the heat arrives); R(t) =
%       sum_k post_full(k) * mu_US(k) from the full posterior (percept
%       after cues and heat are combined).
%    4. MAP assignment z(t) = argmax of the full posterior; sufficient
%       stats of that cause updated.
%
%  Pseudocount a = 0.1 (default): with a = 1 an empty cause is as sharp as
%  a cause holding one observation, centered at mu0 on every feature; the
%  morph (US 0.5, cues 0.5/0.5) matches that template exactly and never
%  reaches the conditioning causes. With a << 1 the new-cause predictive
%  variance is sigma^2 * (1 + 1/a), i.e. broad.

    if nargin < 2; opts = struct(); end
    opts = set_defaults(opts);
    results.opts = opts;

    [T, D] = size(X);
    K = opts.K;
    if opts.alpha == 0; K = 1; end

    Nk = zeros(1, K); SumF = zeros(K, D); z_prev = 1;

    sigma2      = zeros(1, D);
    sigma2(1)   = opts.sigma_us^2;   % SD -> variance
    sigma2(2:D) = opts.sigma_cs^2;

    mu0 = opts.mu0;
    a   = opts.a;

    post = zeros(T, K);
    z    = zeros(T, 1);
    V    = zeros(T, 1);
    R    = zeros(T, 1);

    for t = 1:T
        xt = X(t, :);

        if opts.alpha > 0
            prior = Nk;
            prior(z_prev) = prior(z_prev) + opts.stickiness;
            idx_new = find(prior == 0, 1);
            if isempty(idx_new); idx_new = K; end   % capacity full: reuse last
            prior(idx_new) = opts.alpha;
            prior = prior / sum(prior);

            lik = zeros(1, K); lik_cs = zeros(1, K);
            for k = 1:K
                ll_full = 0; ll_cs = 0;
                for d = 1:D
                    if Nk(k) > 0
                        mu_hat = (SumF(k,d) + a*mu0) / (Nk(k) + a);
                    else
                        mu_hat = mu0;
                    end
                    var_pred = sigma2(d) * (1 + 1/(Nk(k) + a));
                    ll_d = -0.5*log(2*pi*var_pred) - 0.5*(xt(d)-mu_hat)^2/var_pred;
                    ll_full = ll_full + ll_d;
                    if d >= 2; ll_cs = ll_cs + ll_d; end
                end
                lik(k)    = exp(ll_full);
                lik_cs(k) = exp(ll_cs);
            end
            pc = prior .* lik_cs; post_cs   = pc / max(sum(pc), eps);
            pf = prior .* lik;    post_full = pf / max(sum(pf), eps);
        else
            post_cs   = zeros(1,K); post_cs(1)   = 1;
            post_full = zeros(1,K); post_full(1) = 1;
        end

        % read-outs: V from the CS-only posterior (pre-US), R from the
        % full posterior (after cues and heat are combined)
        V_t = 0; R_t = 0;
        for k = 1:K
            if Nk(k) > 0
                mu_us_k = (SumF(k,1) + a*mu0) / (Nk(k) + a);
            else
                mu_us_k = mu0;
            end
            V_t = V_t + post_cs(k)   * mu_us_k;
            R_t = R_t + post_full(k) * mu_us_k;
        end
        V(t) = V_t;
        R(t) = R_t;
        post(t,:) = post_full;

        % MAP assignment + sufficient-stat update
        [~, z_map] = max(post_full);
        Nk(z_map)     = Nk(z_map) + 1;
        SumF(z_map,:) = SumF(z_map,:) + xt;
        z(t) = z_map;  z_prev = z_map;
    end

    results.V    = V;
    results.R    = R;
    results.post = post;          % [T x K], columns = cause ids (NOT reindexed)
    results.z    = z;             % values in 1..K, consistent with post columns
    results.Nk   = Nk;
    results.SumF = SumF;
    results.K    = nnz(Nk > 0);
    results.mu0  = mu0;
    results.a    = a;

end

function opts = set_defaults(opts)
    d = struct('alpha',1,'stickiness',0,'lambda',0,'a',0.1,'mu0',0.5, ...
               'K',10,'sigma_cs',0.2,'sigma_us',0.1);
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(opts,f{i}) || isempty(opts.(f{i})); opts.(f{i}) = d.(f{i}); end
    end
end
