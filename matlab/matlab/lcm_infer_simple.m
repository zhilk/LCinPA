function results = lcm_infer_simple(X, sus, scs, c0, a)
%LCM_INFER_SIMPLE  Local MAP inference for the simple three-cause model,
%  one block of trials.
%
%  Identical inference loop to lcm_infer_percept.m except for the prior:
%  K = 3 fixed slots, prior proportional to visit counts Nk with every
%  empty slot at a fixed weight c0 (no CRP, no alpha). Each trial:
%    1. prior over slots from counts (empty slots at c0), normalized;
%    2. Gaussian feature likelihood of x_t = [US, x_face, x_house] under
%       each slot, means (SumF + a*mu0)/(Nk + a), variances
%       sigma^2 * (1 + 1/(Nk + a)) with sigma = [sus, scs, scs];
%    3. full posterior; read-out R(t) = sum_k post(k) * mu_US(k);
%    4. MAP assignment z(t), sufficient-stat update of that slot.
%
%  Defaults: scs = 0.2, c0 = 0.1, a = 0.1 (fixed in the paper; sweeps in
%  results/). Returns R, post [T x 3], z, Nk, SumF.

    if nargin < 3, scs = 0.2; end
    if nargin < 4, c0  = 0.1; end
    if nargin < 5, a   = 0.1; end

    [T, D] = size(X);
    K = 3; mu0 = 0.5;
    sigma2      = zeros(1, D);
    sigma2(1)   = sus^2;             % SD -> variance
    sigma2(2:D) = scs^2;
    Nk = zeros(1,K); SumF = zeros(K,D);
    R = zeros(T,1); z = zeros(T,1); post = zeros(T,K);
    for t = 1:T
        xt = X(t,:);
        prior = Nk;
        prior(Nk == 0) = c0;         % empty slots at fixed weight
        prior = prior / sum(prior);
        lik = zeros(1,K);
        for k = 1:K
            ll_full = 0;
            for d = 1:D
                if Nk(k) > 0
                    mu_hat = (SumF(k,d) + a*mu0) / (Nk(k) + a);
                else
                    mu_hat = mu0;
                end
                var_pred = sigma2(d) * (1 + 1/(Nk(k) + a));
                ll_full = ll_full - 0.5*log(2*pi*var_pred) ...
                                  - 0.5*(xt(d)-mu_hat)^2/var_pred;
            end
            lik(k) = exp(ll_full);
        end
        pf = prior .* lik;
        post_full = pf / max(sum(pf), eps);

        % US estimate from the full posterior
        R_t = 0;
        for k = 1:K
            if Nk(k) > 0
                mu_us_k = (SumF(k,1) + a*mu0) / (Nk(k) + a);
            else
                mu_us_k = mu0;
            end
            R_t = R_t + post_full(k) * mu_us_k;
        end
        R(t) = R_t;
        post(t,:) = post_full;

        % MAP assignment + sufficient-stat update
        [~, z_map] = max(post_full);
        Nk(z_map)     = Nk(z_map) + 1;
        SumF(z_map,:) = SumF(z_map,:) + xt;
        z(t) = z_map;
    end
    results.R = R; results.post = post; results.z = z;
    results.Nk = Nk; results.SumF = SumF;
end
