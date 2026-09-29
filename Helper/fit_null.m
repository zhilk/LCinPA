function results = fit_null(subT)

%  Null model: ratings driven by the delivered stimulus alone.
%      R(t) = US(t)
%  No expectation, no learning, no latent causes.  The affine observation
%  model (b0, b1, sigma) is the only fitted part, so this is the baseline
%  against which any expectation term must earn its parameters.

    subT = sortrows(subT, 'TrialGlobal');
    US   = subT.TargetVAS / 100;
    CR   = subT.VASRating;
    nTrials = height(subT);

    R = US;                                  % prediction = delivered pain
    [LL, b0, b1, sigma] = rescaled_LL(R, CR);

    results.R      = R;
    results.V      = nan(nTrials,1);         % no expectation term
    results.CRpred = b0 + b1 * R;
    results.LL     = LL;
    results.nPar   = 3;                      % b0, b1, sigma
    results.BIC    = -2*LL + results.nPar * log(nTrials);
    results.b0 = b0; results.b1 = b1; results.sigma = sigma;
end