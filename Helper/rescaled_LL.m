function [LL, b0, b1, sigma] = rescaled_LL(pred, CR)
%RESCALED_LL  Gaussian log-likelihood of ratings under an affine observation model.
%
%  Maps a model's raw prediction onto the observed VAS scale via OLS,
%      CRpred = b0 + b1 * pred,
%  and returns the log-likelihood of CR under Gaussian residuals.
%  b0 absorbs individual rating offset (calibration), b1 the scale gain,
%  sigma is the ML residual SD.
%
%  INPUTS:
%    pred – [n x 1] model prediction (any scale)
%    CR   – [n x 1] observed ratings
%
%  OUTPUTS:
%    LL    – log-likelihood
%    b0,b1 – OLS intercept and slope
%    sigma – ML estimate of residual SD

    n  = numel(CR);
    Xr = [ones(n,1), pred(:)];
    b  = Xr \ CR(:);
    b0 = b(1);
    b1 = b(2);
    resid = CR(:) - Xr*b;
    sigma = sqrt(mean(resid.^2));
    LL    = -0.5*n*log(2*pi) - n*log(sigma) - 0.5*sum((resid/sigma).^2);
end