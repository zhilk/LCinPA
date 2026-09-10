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

% fit random seed for reproducibility
% add particle filter back in
% set a high number of particles to make sure PF likelihood is smooth


    if nargin < 2; opts = struct(); end
    opts = set_defaults(opts);
    results.opts = opts;

    [T, D] = size(X);
    K = opts.K;
    M = opts.M;
    if opts.alpha == 0; K = 1; end

    % Sufficient statistics per cause
    Nk    = zeros(M, K);          % observation count per cause
    SumF  = zeros(M, K, D);          % running sum of features per cause
    % SumF2 = zeros(M, K, D);          % running sum of squared features
    z_prev = ones(M, 1);          % 

    % Feature-specific variance
    sigma2      = ones(1, D);
    sigma2(1)   = opts.sigma_us;  % US variance
    sigma2(2:D) = opts.sigma_cs;  % CS variance

    post_cs = zeros(M,K); 
    post_full = zeros(M,K);

    % Prior mean for all features
    mu0 = 0.5;
    a_pseudo = opts.a;            % pseudocount for prior

    % Storage
    results.V    = zeros(T, 1);
    results.post = zeros(T, K);
    results.z    = zeros(T, 1);
    % z_prev = 1;                   % most recently active cause

    for t = 1:T
        xt = X(t, :);             % [1 x D]: US, CS1, CS2
         
        % % --- Temporal decay of cause counts 
        % if opts.lambda > 0
        %     Nk = exp(-opts.lambda) * Nk;
        % end

        if opts.alpha > 0

            for m = 1:M

                % --- CRP prior  for particle m ---
                prior = Nk(m, :);
                prior(m, z_prev(m)) = prior(z_prev(m)) + opts.stickiness;
                idx_new = find(prior == 0, 1);
    
                if isempty(idx_new) % Expand capacity
                    idx_new = K;
                    % K = K + 1;
                    % Nk(:,K)=0; SumF(:,K,:)=0; SumF2(:,K,:)=0;
                    % results.post(:,K) = 0;
                    % prior(K) = 0; idx_new = K;
                end
                prior(idx_new) = opts.alpha;
                prior = prior / sum(prior);
    
                % --- Likelihood for particle ---
                lik = zeros(1, K); lik_cs = zeros(1, K);
                for k = 1:K
                    log_lik_full = 0; log_lik_cs   = 0;
                    for d = 1:D
                        if Nk(m, k) > 0
                            mu_hat = (SumF(m, k,d) + a_pseudo * mu0) / (Nk(m, k) + a_pseudo);
                        else
                            mu_hat = mu0;
                        end
                        var_pred = sigma2(d) + sigma2(d) / (Nk(m,k) + a_pseudo);
                        ll_d = -0.5 * log(2*pi*var_pred) - 0.5 * (xt(d) - mu_hat)^2 / var_pred;
                        log_lik_full = log_lik_full + ll_d;
                        if d >= 2, log_lik_cs = log_lik_cs + ll_d; end
                    end
                    lik(k)    = exp(log_lik_full);
                    lik_cs(k) = exp(log_lik_cs);
                end

            % --- CS-only posterior (for prediction) ---
            pc = prior .* lik_cs;
            post_cs(m,:) = pc / max(sum(pc), eps);

            % --- Full posterior (CS + US, for cause assignment) ---
            pf = prior .* lik;
            post_full(m,:) = pf / max(sum(pf), eps);
            
            end 
        else
            % alpha = 0: single cause, no inference needed
            post_cs   = zeros(M, K);  post_cs(:,1)   = 1;
            post_full = zeros(M, K);  post_full(:,1)  = 1;
        end

        % --- US prediction: posterior-weighted average of cause means ---
        V_t = 0;
        for m = 1:M 
            for k = 1:K
                if Nk(m, k) > 0
                    mu_us_k = (SumF(m,k,1) + a_pseudo * mu0) / (Nk(m, k) + a_pseudo);
                else
                    mu_us_k = mu0;
                end
                V_t = V_t + post_cs(m,k) * mu_us_k;
            end
        end 
        results.V(t) = V_t/M; %averaged over particles
        results.post(t, 1:K) = post_full;

        % --- MAP assignment ---
        if M==1 
        [~, z_map] = max(post_full(1,:));
        Nk(1, z_map)       = Nk(1, z_map) + 1;
        SumF(1, z_map, :)  = SumF(1, z_map, :) + reshape(xt,1,1,D);
        results.z(t) = z_map;
        z_prev(1) = z_map;
        else
            NkOld = Nk; SumFOld = SumF;
            pm = sum(post_full, 2);                    % M×1
            for m = 1:M
                row = find(rand()*sum(pm) < cumsum(pm), 1);
                Nk(m,:)     = NkOld(row,:);
                SumF(m,:,:) = SumFOld(row,:,:);
                pr  = post_full(row,:) / sum(post_full(row,:));
                col = find(rand() < cumsum(pr), 1);
                Nk(m,col)     = Nk(m,col) + 1;
                SumF(m,col,:) = SumF(m,col,:) + reshape(xt,1,1,D);
                z_prev(m)     = col;
            end
            results.z(t) = mode(z_prev);
        end 
        
    end

    % Trim unused causes
    active = sum(results.post, 1) > 0;
    results.post = results.post(:, active);
    results.K = sum(active);

    % Store final sufficient statistics (for cause assignment analysis)
    results.Nk   = mean(Nk, 1);
    results.Nk   = results.Nk(active);
    results.SumF = squeeze(mean(SumF, 1));      % M×K×D -> K×D
    results.SumF = results.SumF(active, :);
    results.mu0  = mu0;
    results.a_pseudo = a_pseudo;
end


%% ========================================================================
function opts = set_defaults(opts)
    if ~isfield(opts,'alpha')      || isempty(opts.alpha);      opts.alpha = 1;     end
    if ~isfield(opts,'stickiness') || isempty(opts.stickiness); opts.stickiness = 0; end
    if ~isfield(opts,'lambda')     || isempty(opts.lambda);     opts.lambda = 0;    end
    if ~isfield(opts,'M')          || isempty(opts.M);          opts.M = 1;     end
    if ~isfield(opts,'a')          || isempty(opts.a);          opts.a = 1;         end
    if ~isfield(opts,'b')          || isempty(opts.b);          opts.b = 1;         end
    if ~isfield(opts,'K')          || isempty(opts.K);          opts.K = 10;        end
    if ~isfield(opts,'sigma_cs')   || isempty(opts.sigma_cs);   opts.sigma_cs = 0.5; end
    if ~isfield(opts,'sigma_us')   || isempty(opts.sigma_us);   opts.sigma_us = 0.5; end
end
