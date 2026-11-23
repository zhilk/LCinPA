function V = rw_forward(CS, US, eta, blocks)
%RW_FORWARD  Rescorla-Wagner cue expectation with per-block naive reset.
%
%      V(t)   = w' * CS(t)
%      w_{t+1} = w + eta * CS(t) * (US(t) - V(t))
%
%  Weights start naive at [0.5 0.5] and reset to naive at every block
%  boundary. Returns the expectation V. For the RW model R(t) = V(t);
%  for RW-BCC, predict_model.m blends it with US: R = kappa*V + (1-kappa)*US.

    n  = numel(US);
    V = zeros(n, 1);
    w  = [0.5; 0.5];
    for t = 1:n
        if t > 1 && blocks(t) ~= blocks(t-1)
            w = [0.5; 0.5];                    % naive reset at block boundary
        end
        x = CS(t, :)';
        V(t) = w' * x;
        w = w + eta * x * (US(t) - V(t));
    end
end
