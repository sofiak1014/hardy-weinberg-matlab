function [chi2, pval, pHat, Fhat] = hw_chisq(counts)
% HW_CHISQ  Chi-square goodness-of-fit test for Hardy-Weinberg equilibrium.
%
%   [chi2, pval, pHat, Fhat] = hw_chisq([nAA nAa naa])
%
%   Steps:
%     1) Estimate p from the data:  pHat = (2*nAA + nAa) / (2n)
%     2) Expected counts under HWE: n * [p^2, 2pq, q^2]
%     3) chi2 = sum( (observed - expected)^2 / expected )
%     4) p-value from a chi-square distribution with 1 degree of freedom
%        (3 genotype classes - 1 - 1 because p was estimated from the same data)
%
%   Fhat = 1 - (observed Aa) / (expected Aa); positive means too few heterozygotes.
%   No continuity correction is applied.
%   Returns NaN if the sample is monomorphic (only one allele present).

    n    = sum(counts);
    pHat = (2*counts(1) + counts(2)) / (2*n);
    qHat = 1 - pHat;

    if pHat == 0 || pHat == 1
        chi2 = NaN; pval = NaN; Fhat = NaN;
        return
    end

    expected = n * [pHat^2, 2*pHat*qHat, qHat^2];
    chi2     = sum((counts - expected).^2 ./ expected);

    % Upper tail of chi-square with 1 df:  P(X > x) = gammainc(x/2, 1/2, 'upper')
    % (base MATLAB; chi2cdf would need the Statistics Toolbox)
    pval = gammainc(chi2/2, 0.5, 'upper');
    Fhat = 1 - counts(2) / expected(2);
end
