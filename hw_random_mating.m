function [counts, pParents, pOffspring] = hw_random_mating(genoFreq, N)
% HW_RANDOM_MATING  One generation of random mating.
%
%   [counts, pParents, pOffspring] = hw_random_mating(genoFreq, N)
%
%   genoFreq  parent genotype frequencies [fAA fAa faa] (must sum to 1)
%   N         number of offspring to produce
%
%   counts      offspring genotype counts [nAA nAa naa]
%   pParents    frequency of allele A among the parents
%   pOffspring  frequency of allele A among the offspring
%
% How it works (step by step):
%   1) Parents' gene pool: p = fAA + fAa/2   (every AA gives two A's, every Aa gives one)
%   2) Each offspring gets TWO gametes, each one is A with probability p, chosen independently
%      (that independence IS random mating)
%   3) Count how many A's each offspring got: 2 = AA, 1 = Aa, 0 = aa

    if abs(sum(genoFreq) - 1) > 1e-9
        error('genoFreq must sum to 1.');
    end
    pParents = genoFreq(1) + genoFreq(2)/2;

    gametes  = rand(N, 2) < pParents;      % N x 2 matrix of true/false (true = A)
    nA       = sum(gametes, 2);            % number of A alleles per offspring (0, 1, or 2)

    counts     = [sum(nA == 2), sum(nA == 1), sum(nA == 0)];
    pOffspring = (2*counts(1) + counts(2)) / (2*N);
end
