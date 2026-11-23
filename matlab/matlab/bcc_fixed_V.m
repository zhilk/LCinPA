function V = bcc_fixed_V(sub)
%BCC_FIXED_V  Cue expectation with fixed objective weights (pure Bayesian
%  cue combination, no learning).
%
%      V(t) = CS(t) * w_obj(block)
%
%  w_obj extracts CS->US mapping from from the first conditioning trial of 
%  that cue. V is therefore objectively correct cue-based expectation from 
%  trial one as CS->US mapping is constant in conditioning.
%
%  predict_model.m blends it with US: R = kappa*V + (1-kappa)*US for bcc fit.


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
