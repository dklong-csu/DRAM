function C = adaptCov(dim, samps, eps)
    sd = 2.4^2/dim;
    CovN = cov(samps);
    C = sd*CovN + sd*eps*eye(dim);
end