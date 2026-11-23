function g = model_grids()
%MODEL_GRIDS  All parameter grids and free-parameter counts in one place.
%
%  eta     RW learning rate, 50 values (fit_rw)
%  eta30   RW learning rate, 30 values (fit_rwbcc joint grid)
%  kappa   expectation weight, 30 values (fit_bcc, fit_rwbcc)
%  sus     sigma_us, 12 values from near-deterministic pain routing (0.005)
%          to pain nearly uninformative for cause assignment (0.40);
%          shared by the simple and the percept model
%  alpha   CRP concentration, 12 values; 0 = single cause, values above
%          ~10 are behaviorally saturated at these sequence lengths
%  npar    free parameters per model, including b0, b1, sigma

    g.eta   = linspace(0.01, 1, 50);
    g.eta30 = linspace(0.01, 1, 30);
    g.kappa = linspace(0, 1, 30);
    g.sus   = [0.005 0.01 0.015 0.02 0.03 0.05 0.08 0.10 0.12 0.16 0.25 0.40];
    g.alpha = [0 0.1 0.3 1 2 3 5 8 10 15 25 50];
    g.npar  = struct('null',3, 'rw',4, 'bcc',4, 'rwbcc',5, ...
                     'simple',4, 'percept',5);
end
