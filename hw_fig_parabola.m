% HW_FIG_PARABOLA  Figure 2: the Hardy-Weinberg curves, with simulated data on top.
%
% Left panel : genotype frequencies (AA, Aa, aa) vs allele frequency p.
%              Lines = theory (p^2, 2pq, q^2). Dots = simulated samples.
%              Notice Aa peaks at 0.5 when p = 0.5.
% Right panel: heterozygote (Aa) frequency vs p for different inbreeding
%              coefficients F. Higher F pulls the curve down (fewer heterozygotes).
%
% Needs hw_expected.m and hw_sample.m in the same folder.

clear; clc; close all; rng(3);

n      = 1000;                    % individuals per simulated sample
pGrid  = linspace(0, 1, 201);     % smooth grid for the theory lines
pPts   = 0.05:0.05:0.95;          % allele frequencies where we simulate
col    = lines(3);                % three plot colors

figure('Position', [100 100 1000 430]);

% ---------- Left panel: Hardy-Weinberg curves ----------
subplot(1, 2, 1); hold on;
theory = zeros(numel(pGrid), 3);
for i = 1:numel(pGrid)
    theory(i, :) = hw_expected(pGrid(i));
end
obs = zeros(numel(pPts), 3);
for i = 1:numel(pPts)
    obs(i, :) = hw_sample(n, pPts(i)) / n;     % counts -> frequencies
end
names = {'AA', 'Aa', 'aa'};
for k = 1:3
    plot(pGrid, theory(:, k), '-', 'Color', col(k, :), 'LineWidth', 1.8, ...
         'DisplayName', [names{k} ' (theory)']);
end
for k = 1:3
    plot(pPts, obs(:, k), 'o', 'Color', col(k, :), 'MarkerFaceColor', col(k, :), ...
         'MarkerSize', 5, 'HandleVisibility', 'off');
end
plot(nan, nan, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, ...
     'DisplayName', sprintf('simulated (n = %d)', n));
hold off;
xlabel('Allele frequency p (frequency of A)');
ylabel('Genotype frequency');
title('Hardy-Weinberg: genotype frequencies vs p');
legend('Location', 'north');
grid on; xlim([0 1]); ylim([0 1]);

% ---------- Right panel: effect of inbreeding on heterozygotes ----------
subplot(1, 2, 2); hold on;
Fvals = [0 0.3 0.6];
for j = 1:numel(Fvals)
    F = Fvals(j);
    het = zeros(numel(pGrid), 1);
    for i = 1:numel(pGrid)
        f = hw_expected(pGrid(i), F);
        het(i) = f(2);
    end
    plot(pGrid, het, '-', 'Color', col(j, :), 'LineWidth', 1.8, ...
         'DisplayName', sprintf('F = %.1f (theory)', F));
    hetObs = zeros(numel(pPts), 1);
    for i = 1:numel(pPts)
        c = hw_sample(n, pPts(i), F);
        hetObs(i) = c(2) / n;
    end
    plot(pPts, hetObs, 'o', 'Color', col(j, :), 'MarkerFaceColor', col(j, :), ...
         'MarkerSize', 5, 'HandleVisibility', 'off');
end
hold off;
xlabel('Allele frequency p (frequency of A)');
ylabel('Heterozygote (Aa) frequency');
title('Inbreeding lowers heterozygote frequency');
legend('Location', 'north');
grid on; xlim([0 1]); ylim([0 0.6]);

saveas(gcf, 'fig2_parabola.png');
