% HW_FIG_CONVERGENCE  Figure 1: genotype frequencies reach Hardy-Weinberg
% proportions after ONE generation of random mating, then stay there.
%
% Start: 60% AA, 0% Aa, 40% aa  (very far from equilibrium, no heterozygotes).
% Needs hw_random_mating.m and hw_expected.m in the same folder.

clear; clc; close all; rng(2);

start = [0.6 0 0.4];     % starting genotype frequencies [AA Aa aa]
N     = 1e5;             % offspring per generation (large, so drift is tiny)
gens  = 4;               % number of generations of random mating

% Run the generations. Row g holds the genotype frequencies of generation g-1.
G = zeros(gens + 1, 3);
G(1, :) = start;
for g = 1:gens
    counts = hw_random_mating(G(g, :), N);
    G(g + 1, :) = counts / N;          % counts -> frequencies
end

% Hardy-Weinberg prediction from the STARTING allele frequency
p        = start(1) + start(2)/2;      % = 0.6
expected = hw_expected(p);             % [p^2 2pq q^2] = [0.36 0.48 0.16]

% Print a table so you can read the numbers
pAllele = G(:, 1) + G(:, 2)/2;
fprintf('Gen    AA       Aa       aa      p(A)\n');
for g = 0:gens
    fprintf('%2d   %.4f   %.4f   %.4f   %.4f\n', g, G(g+1, :), pAllele(g+1));
end
fprintf('HW   %.4f   %.4f   %.4f   %.4f  (expected)\n', expected, p);

% Plot
figure('Position', [100 100 760 460]);
bar(0:gens, G, 'grouped');
hold on;
for k = 1:3     % dashed lines = Hardy-Weinberg expectation for AA, Aa, aa
    plot([-0.5 gens + 0.5], [expected(k) expected(k)], 'k--', ...
         'HandleVisibility', 'off');
end
plot(nan, nan, 'k--');   % one legend entry for the dashed lines
hold off;
xlabel('Generation');
ylabel('Genotype frequency');
title('One generation of random mating reaches Hardy-Weinberg proportions');
legend({'AA', 'Aa', 'aa', 'HW expectation (p^2, 2pq, q^2)'}, 'Location', 'north');
ylim([0 0.75]);
grid on;

saveas(gcf, 'fig1_convergence.png');
