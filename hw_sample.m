function counts = hw_sample(n, p, F)
% HW_SAMPLE  Draw a random sample of n individuals and count genotypes.
%
%   counts = hw_sample(n, p)      sample from a Hardy-Weinberg population
%   counts = hw_sample(n, p, F)   sample from a population with inbreeding coefficient F
%
%   counts = [nAA nAa naa]. Uses base MATLAB only (no toolboxes).

    if nargin < 3, F = 0; end
    f = hw_expected(p, F);
    u = rand(n, 1);                         % one random number per individual
    counts = [sum(u < f(1)), ...
              sum(u >= f(1) & u < f(1) + f(2)), ...
              sum(u >= f(1) + f(2))];
end
