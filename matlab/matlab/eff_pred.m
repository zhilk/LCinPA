function e = eff_pred(sub, R, b0, b1)
%EFF_PRED  Model-predicted morph effect in VAS points.
%
%  Mean predicted rating (b0 + b1*R) on morph trials, nocebo block minus
%  placebo block. Positive = predicted context effect in the observed
%  direction (morphs rated higher in the nocebo block).

    pr = b0 + b1*R;
    e = mean(pr(sub.isMix & sub.blocks == "nocebo")) ...
      - mean(pr(sub.isMix & sub.blocks == "placebo"));
end
