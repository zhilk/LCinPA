function f = fit_on_grid(pred, params, CR, mask, nPar)
%FIT_ON_GRID  Grid maximum likelihood of the rating model over all columns
%  of a prediction matrix.
%
%  For every grid column R the rating model CR = b0 + b1*R + N(0, sigma^2)
%  is solved in closed form on the scored trials (mask): b0, b1 by OLS,
%  sigma as the residual RMS; then LL = -n/2*log(2*pi) - n*log(sigma) - n/2
%  (the quadratic term equals n/2 at the ML sigma). Identical to
%  rescaled_LL.m applied per column, vectorized. Columns with (numerically)
%  constant predictions get b1 = 0, so they reduce to an intercept-only
%  fit. Ties resolve to the first column (see predict_model.m).
%
%  Returns struct f: idx (winning column), p (its parameter row), LL, b0,
%  b1, sigma, BIC with n = nnz(mask).

    P = pred(mask, :);
    y = CR(mask);
    n = numel(y);

    mx  = mean(P, 1);
    my  = mean(y);
    vx  = mean(P.^2, 1) - mx.^2;
    cxy = (y' * P)/n - mx*my;
    b1  = cxy ./ max(vx, 1e-12);
    b1(vx <= 1e-12) = 0;
    b0  = my - b1 .* mx;
    res = y - (b0 + P .* b1);
    sig = max(sqrt(mean(res.^2, 1)), 1e-6);
    LL  = -0.5*n*log(2*pi) - n*log(sig) - 0.5*n;

    [f.LL, f.idx] = max(LL);
    f.p     = params(f.idx, :);
    f.b0    = b0(f.idx);
    f.b1    = b1(f.idx);
    f.sigma = sig(f.idx);
    f.BIC   = -2*f.LL + nPar*log(n);
end
