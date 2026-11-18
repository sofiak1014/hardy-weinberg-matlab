function f = hw_expected(p, F)
% HW_EXPECTED  Expected genotype frequencies [AA Aa aa] for allele freq p.
%
%   f = hw_expected(p)      random mating (Hardy-Weinberg): [p^2, 2pq, q^2]
%   f = hw_expected(p, F)   with inbreeding coefficient F (F = 0 is random mating)
%
%   AA = p^2 + F*p*q,   Aa = 2*p*q*(1 - F),   aa = q^2 + F*p*q
%   (inbreeding form: Weir & Goudet, Summer Institute in Statistical Genetics)

    if nargin < 2, F = 0; end
    q = 1 - p;
    f = [p^2 + F*p*q,  2*p*q*(1 - F),  q^2 + F*p*q];
end
