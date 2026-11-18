% HW_TESTS  Checks that the Hardy-Weinberg functions behave correctly.
% Run this file. Every line should print PASS.
clear; clc; rng(1);

msgs  = {'*** FAIL ***', 'PASS'};
check = @(name, ok) fprintf('%-62s %s\n', name, msgs{logical(ok) + 1});

% 1) Formula checks (done by hand)
f = hw_expected(0.5);
check('p = 0.5, F = 0  ->  [0.25 0.5 0.25]', all(abs(f - [0.25 0.5 0.25]) < 1e-12));
check('genotype frequencies sum to 1 (p = 0.3, F = 0.2)', abs(sum(hw_expected(0.3, 0.2)) - 1) < 1e-12);
f = hw_expected(0.5, 1);
check('F = 1 (full inbreeding): no heterozygotes', abs(f(2)) < 1e-12);

% 2) Chi-square on hand-calculated cases
%    [40 0 60]: p = 0.4, expected = [16 48 36], chi2 = 36 + 48 + 16 = 100
[c, pv] = hw_chisq([40 0 60]);
check('counts [40 0 60]: chi2 = 100 (by hand)', abs(c - 100) < 1e-9);
check('counts [40 0 60]: p-value about 1.52e-23', abs(pv/1.5239706e-23 - 1) < 1e-3);
[c, pv] = hw_chisq([4 0 6]);
check('counts [4 0 6]: chi2 = 10 (same proportions, n = 10)', abs(c - 10) < 1e-9);
check('counts [4 0 6]: p-value about 0.001565', abs(pv - 0.0015654) < 1e-6);
[c, pv, pHat, Fhat] = hw_chisq([25 50 25]);
check('perfect HWE counts [25 50 25]: chi2 = 0, p-value = 1, F = 0', ...
      abs(c) < 1e-12 && abs(pv - 1) < 1e-12 && abs(Fhat) < 1e-12);

% 3) One generation of random mating (start with NO heterozygotes)
[counts, pPar, pOff] = hw_random_mating([0.6 0 0.4], 1e6);
obs = counts / sum(counts);
check('allele frequency unchanged after 1 generation', abs(pOff - pPar) < 0.005);
check('offspring genotypes reach [p^2 2pq q^2] = [0.36 0.48 0.16]', ...
      all(abs(obs - [0.36 0.48 0.16]) < 0.005));

% 4) False-positive rate of the test on TRUE Hardy-Weinberg data
%    (expect about 5% at alpha = 0.05 for a common allele and a decent sample size)
reps = 5000; rej = 0;
for r = 1:reps
    [~, pv] = hw_chisq(hw_sample(500, 0.5));
    rej = rej + (pv < 0.05);
end
rate = rej / reps;
fprintf('   false-positive rate (p = 0.5, n = 500): %.4f\n', rate);
check('false-positive rate within 0.04 to 0.06', rate > 0.04 && rate < 0.06);

% 5) The test should DETECT inbreeding
reps = 1000; rej = 0;
for r = 1:reps
    [~, pv] = hw_chisq(hw_sample(500, 0.5, 0.2));
    rej = rej + (pv < 0.05);
end
power = rej / reps;
fprintf('   detection rate (F = 0.2, n = 500): %.4f\n', power);
check('detects F = 0.2 in at least 95% of samples', power > 0.95);
