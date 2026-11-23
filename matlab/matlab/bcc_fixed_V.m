function V = bcc_fixed_V(sub)
%BCC_FIXED_V  Cue expectation with fixed objective weights (pure Bayesian
%  cue combination, no learning).
%
%      V(t) = CS(t) * w_obj(block)
%
%  w_obj holds the delivered heat per cue per block (read from the first
%  conditioning trial of each pure-cue category; delivery is constant
%  within a block's conditioning phase), so V is the normatively correct
%  cue-based expectation from trial one. fit_bcc blends it with US:
%  R = kappa*V + (1-kappa)*US.

    V = zeros(sub.n, 1);
    for bi = 1:numel(sub.ub)
        b = sub.ub(bi);
        m = sub.blocks == b;
        w = zeros(2, 1);
        w(1) = sub.US(find(sub.isCond & m & sub.cat == "face",  1));
        w(2) = sub.US(find(sub.isCond & m & sub.cat == "house", 1));
        V(m) = sub.CS(m, :) * w;
    end
end
