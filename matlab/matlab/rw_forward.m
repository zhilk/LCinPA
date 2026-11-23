function mu = rw_forward(CS, US, eta, blocks)
%RW_FORWARD  Rescorla-Wagner cue expectation with per-block naive reset.
%
%      mu(t)   = w' * CS(t)
%      w_{t+1} = w + eta * CS(t) * (US(t) - mu(t))
%
%  Weights start naive at [0.5 0.5] and reset to naive at every block
%  boundary (the contingencies are re-conditioned per block; all models in
%  the set begin each block naive). Returns the expectation mu, i.e. the
%  pure RW read-out; fit_rwbcc blends it with US downstream.

    n  = numel(US);
    mu = zeros(n, 1);
    w  = [0.5; 0.5];
    for t = 1:n
        if t > 1 && blocks(t) ~= blocks(t-1)
            w = [0.5; 0.5];                    % naive reset at block boundary
        end
        x = CS(t, :)';
        mu(t) = w' * x;
        w = w + eta * x * (US(t) - mu(t));
    end
end
