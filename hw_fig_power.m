% HW_FIG_POWER  Figure 3: how well does the chi-square HWE test work?
%
% Left panel : DETECTION RATE vs sample size for different inbreeding levels F
%              (p = 0.5). F = 0 is a true Hardy-Weinberg population, so that
%              line is the false-positive rate and should sit near 0.05.
% Right panel: FALSE-POSITIVE RATE vs sample size for common vs rare alleles
%              (true HWE, F = 0). Should be about 0.05 everywhere if the test
%              behaved perfectly. Wigginton, Cutler & Abecasis (2005) warn the
%              chi-square test can misbehave when an allele is rare.
%
% Needs hw_sample.m, hw_expected.m, hw_chisq.m in the same folder.
% Run time: roughly a minute. Lower 'reps' to go faster (results get noisier).

clear; clc; close all; rng(4);

reps  = 1000;                          % simulated samples per point
alpha = 0.05;                          % significance level
nVals = [25 50 100 200 400 800];       % sample sizes

col = lines(4);
figure('Position', [100 100 1000 430]);

% ---------- Left: detection rate vs n, for several F (p = 0.5) ----------
subplot(1, 2, 1); hold on;
Fvals = [0 0.1 0.2 0.3];
for j = 1:numel(Fvals)
    rate = zeros(size(nVals));
    for i = 1:numel(nVals)
        rej = 0;
        for r = 1:reps
            [~, pv] = hw_chisq(hw_sample(nVals(i), 0.5, Fvals(j)));
            rej = rej + (pv < alpha);      % NaN (monomorphic) counts as "not rejected"
        end
        rate(i) = rej / reps;
    end
    plot(nVals, rate, '-o', 'Color', col(j, :), 'LineWidth', 1.6, ...
         'MarkerFaceColor', col(j, :), 'DisplayName', sprintf('F = %.1f', Fvals(j)));
end
plot([nVals(1) nVals(end)], [alpha alpha], 'k:', 'HandleVisibility', 'off');
hold off;
set(gca, 'XScale', 'log'); xticks(nVals); xlim([nVals(1) nVals(end)]);
xlabel('Sample size n (individuals)');
ylabel('Fraction of samples flagged (p-value < 0.05)');
title('Detecting inbreeding (p = 0.5)');
legend('Location', 'southeast'); grid on; ylim([0 1]);

% ---------- Right: false-positive rate vs n, common vs rare allele ----------
subplot(1, 2, 2); hold on;
pVals = [0.5 0.1 0.02];
for j = 1:numel(pVals)
    rate = zeros(size(nVals));
    for i = 1:numel(nVals)
        rej = 0;
        for r = 1:reps
            [~, pv] = hw_chisq(hw_sample(nVals(i), pVals(j), 0));
            rej = rej + (pv < alpha);
        end
        rate(i) = rej / reps;
    end
    plot(nVals, rate, '-o', 'Color', col(j, :), 'LineWidth', 1.6, ...
         'MarkerFaceColor', col(j, :), 'DisplayName', sprintf('p = %.2f', pVals(j)));
end
plot([nVals(1) nVals(end)], [alpha alpha], 'k--', 'DisplayName', 'target (0.05)');
hold off;
set(gca, 'XScale', 'log'); xticks(nVals); xlim([nVals(1) nVals(end)]);
xlabel('Sample size n (individuals)');
ylabel('False-positive rate (true HWE)');
title('False positives: common vs rare allele');
legend('Location', 'north'); grid on; ylim([0 0.15]);

saveas(gcf, 'fig3_power.png');
