%% Figure 6: Economically optimal upgrading and carbon-account safety
% Standalone MATLAB script: no user-defined or local functions.
%
% Scientific purpose
% ------------------
% Panel (a) reports the preferred major-upgrade regime WITHIN a finite
% candidate-policy library on a conditional R-A state slice.  This is not
% the global optimum of the unrestricted hybrid-control problem.
%
% Panel (b) identifies the non-mechanical deviation region in which the
% current balance is nonnegative and an upgrade is selected within the
% policy library, but the remaining-horizon account path becomes negative.
%
% Numerical interpretation
% ------------------------
% This is a fast policy-library approximation for regional visualization,
% not a replacement for the continuous-time first-order conditions in the
% theoretical analysis.  It avoids repeated fmincon calls at every state
% point and normally runs within seconds.
%
% The figure uses a short-horizon decision slice corresponding to t = 8
% and T-t = 2.  kappaSlice is a conditional current state, rather than a
% change in the structural model parameters.

clearvars;
clc;
close all;

tic;

%% 1. Structural parameters: consistent with the benchmark simulations

Tfull = 10.0;
tSlice = 8.0;
H = Tfull - tSlice;
N = 40;
dt = H / N;

r = 0.05;
y = 1.00;
price = 3.00;
pE = 0.60;
chi = 0.35;
Q = 0.50;
deltaA = 0.02;

kappaBar = 1.00;
kappaSlice = 0.85;

beta = 0.25;
deltaR = 0.05;
etaR = 0.45;

etaK = 0.35;
deltaK = 0.12;

cm = 0.08;
am = 0.45;
em = 0.03;
em2 = 0.10;

aR = 1.20;
ay = 0.45;

% Major-upgrade cost and effectiveness
K0 = 0.05;       % CHECK against manuscript Table 2 (Fig. 5 used 0.08, unused there)
K1 = 0.25;
K2 = 0.70;

xi1 = 0.14;
xi2 = 0.35;

gammaR = 0.65;
gammaK = 0.65;

vA = 0.50;
phiA = 8.00;
vK = 0.80;
phiR = 1.50;

mbar = 1.00;

%% 2. R-A state grid

nR = 61;
nA = 51;

RGrid = linspace(0.05, 1.15, nR);
% Restrict the initial carbon balance to A0 >= 0.  This removes cases that
% are mechanically unsafe at the beginning of the decision window and
% focuses the plot on endogenous crossings of the zero-balance boundary.
AGrid = linspace(0.00, 1.20, nA);
[RR, AA] = meshgrid(RGrid, AGrid);

%% 3. Fast policy library
% Maintenance is represented by constant effort levels over the remaining
% horizon.  At most one candidate major upgrade is considered on this
% short-horizon slice.  Finer sets improve resolution but increase runtime.

mSet = 0:0.10:mbar;
sigmaSet = [0.25, 0.75, 1.25, 1.75];    % time measured from tSlice
omegaSet = 0.10:0.10:1.00;

candM = [];
candSigma = [];
candOmega = [];

for im = 1:numel(mSet)
    mVal = mSet(im);

    % No major upgrade under the selected maintenance level
    candM(end+1,1) = mVal;
    candSigma(end+1,1) = NaN;
    candOmega(end+1,1) = 0;

    % One candidate major upgrade
    for is = 1:numel(sigmaSet)
        for iw = 1:numel(omegaSet)
            candM(end+1,1) = mVal;
            candSigma(end+1,1) = sigmaSet(is);
            candOmega(end+1,1) = omegaSet(iw);
        end
    end
end

nCandidates = numel(candM);

%% 4. Enumerate candidate policies on the complete state grid

bestValue = -Inf(size(RR));
bestM = zeros(size(RR));
bestSigma = NaN(size(RR));
bestOmega = zeros(size(RR));
bestMinA = NaN(size(RR));
bestTerminalA = NaN(size(RR));
bestBAccount = NaN(size(RR));    % Proposition 3: account-adjusted payback
bestUpgradeAplus = NaN(size(RR));  % account immediately AFTER chosen impulse

discA = exp(-deltaA * dt);
discT = exp(-r * H);

for ic = 1:nCandidates

    m = candM(ic);
    sigma = candSigma(ic);
    omega = candOmega(ic);

    R = RR;
    K = kappaSlice * ones(size(RR));
    A = AA;

    J = zeros(size(RR));
    minAPath = A;
    upgradeApplied = false;
    upgradeAplus = NaN(size(RR));

    % Exact transition coefficients for R and K under constant m
    lambdaR = deltaR + etaR * m;
    Rss = beta * y / lambdaR;

    lambdaK = etaK * m + deltaK;
    Kss = etaK * m * kappaBar / lambdaK;

    for it = 1:N

        tau0 = (it - 1) * dt;
        tauMid = tau0 + 0.5 * dt;
        tau1 = it * dt;

        % Midpoint and endpoint states before any impulse at tau1
        Rmid = Rss + (R - Rss) .* exp(-lambdaR * 0.5 * dt);
        Rnext = Rss + (R - Rss) .* exp(-lambdaR * dt);

        Kmid = Kss + (K - Kss) .* exp(-lambdaK * 0.5 * dt);
        Knext = Kss + (K - Kss) .* exp(-lambdaK * dt);
        Kmid = max(Kmid, 1e-8);

        carbonFlow = Q ...
                   - chi * y ./ Kmid ...
                   - em * m ...
                   - 0.5 * em2 * m^2;

        Anext = A .* discA ...
              + carbonFlow .* (1 - discA) / deltaA;

        operatingCost = pE * y ./ Kmid ...
                      + 0.5 * ay * (y ./ Kmid).^2;

        maintenanceCost = cm * m + 0.5 * am * m^2;
        privacyCost = 0.5 * aR * Rmid.^2;

        currentProfit = price * y ...
                      - operatingCost ...
                      - maintenanceCost ...
                      - privacyCost;

        J = J + exp(-r * tauMid) * dt .* currentProfit;

        R = Rnext;
        K = Knext;
        A = Anext;
        minAPath = min(minAPath, A);

        % Apply the candidate upgrade at the end of this interval
        if ~upgradeApplied && isfinite(sigma) && abs(tau1 - sigma) <= 0.5*dt

            R = (1 - gammaR * omega) .* R;
            K = K + gammaK * omega .* (kappaBar - K);
            A = A - xi1 * omega - 0.5 * xi2 * omega^2;
            upgradeAplus = A;

            upgradeCost = K0 + K1 * omega + 0.5 * K2 * omega^2;
            J = J - exp(-r * sigma) * upgradeCost;

            minAPath = min(minAPath, A);
            upgradeApplied = true;
        end
    end

    % Correct terminal carbon-account value:
    % Phi_A(A_T) = vA*A_T - 0.5*phiA*max(-A_T,0)^2.
    % Using vA*max(A_T,0) would incorrectly ignore the linear loss on A_T<0.
    Aminus = max(-A, 0);
    terminalValue = vA * A ...
                  - 0.5 * phiA * Aminus.^2 ...
                  + vK * K ...
                  - 0.5 * phiR * R.^2;

    J = J + discT .* terminalValue;

    % Proposition 3 diagnostic for this fixed maintenance path. Both
    % counterfactual branches start from the same state at upgrade time;
    % one has the impulse and one does not. The carbon flows Q and
    % maintenance-emission terms therefore cancel in their difference.
    BAccount = NaN;
    if isfinite(sigma) && omega > 0
        uSave = (0:dt:(H-sigma))';
        kBefore = Kss + (kappaSlice-Kss)*exp(-lambdaK*sigma);
        kAfter = kBefore + gammaK*omega*(kappaBar-kBefore);
        kNo = Kss + (kBefore-Kss)*exp(-lambdaK*uSave);
        kYes = Kss + (kAfter-Kss)*exp(-lambdaK*uSave);
        avoided = chi*y./kNo - chi*y./kYes;
        BAccount = trapz(uSave,exp(deltaA*uSave).*avoided) ...
                   - (xi1*omega + 0.5*xi2*omega^2);
    end

    improve = J > bestValue;

    bestValue(improve) = J(improve);
    bestM(improve) = m;
    bestSigma(improve) = sigma;
    bestOmega(improve) = omega;
    bestMinA(improve) = minAPath(improve);
    bestTerminalA(improve) = A(improve);
    bestBAccount(improve) = BAccount;
    bestUpgradeAplus(improve) = upgradeAplus(improve);
end

%% 5. Classification

omegaTol = 1e-8;

% 0 = no major upgrade; 1 = partial upgrade; 2 = full upgrade
upgradeRegime = zeros(size(bestOmega));
upgradeRegime(bestOmega > omegaTol & bestOmega < 1 - omegaTol) = 1;
upgradeRegime(bestOmega >= 1 - omegaTol) = 2;

% Keep the EXACT discrete policy-library classification.  Majority filtering
% can visually reassign a grid cell to a different policy while retaining
% the safety state generated by the original optimal candidate.
upgradeRegimePlot = upgradeRegime;

upgradeSelected = bestOmega > omegaTol;
dynamicSafe = bestMinA >= -1e-10;
deviationRegion = (AA >= 0) & upgradeSelected & ~dynamicSafe;

%% 6. Publication-style figure

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 11);

fig = figure('Name', ...
    'Figure 6: Policy-library upgrading and carbon-account safety', ...
    'Color', 'w', 'Position', [80, 160, 1040, 440]);

tl = tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

% ---------------- Panel (a): economic upgrade regimes ------------------
ax1 = nexttile;

contourf(ax1, RGrid, AGrid, upgradeRegimePlot, ...
    [-0.5, 0.5, 1.5, 2.5], 'LineColor', 'none');

colormap(ax1, [0.93 0.94 0.95; ...
               0.65 0.69 0.73; ...
               0.25 0.28 0.31]);

hold(ax1, 'on');
[~, hSafety1] = contour(ax1, RGrid, AGrid, bestMinA, ...
    [0 0], 'k-', 'LineWidth', 2.0);
hold(ax1, 'off');

xlabel(ax1, 'Privacy risk, $R(t_s)$', 'Interpreter', 'latex');
ylabel(ax1, 'Carbon-account balance, $A(t_s^-)$', 'Interpreter', 'latex');
title(ax1, '(a) Preferred upgrade regimes (policy library)');
xlim(ax1, [RGrid(1), RGrid(end)]);
ylim(ax1, [AGrid(1), AGrid(end)]);
box(ax1, 'on');
grid(ax1, 'on');
set(ax1, 'Layer', 'top', 'LineWidth', 1.0);

hNo = patch(ax1, NaN, NaN, [0.93 0.94 0.95], ...
    'EdgeColor', 'k', 'DisplayName', 'No major upgrade');
hPartial = patch(ax1, NaN, NaN, [0.65 0.69 0.73], ...
    'EdgeColor', 'k', 'DisplayName', 'Partial upgrade');
hFull = patch(ax1, NaN, NaN, [0.25 0.28 0.31], ...
    'EdgeColor', 'k', 'DisplayName', 'Full upgrade');
set(hSafety1, 'DisplayName', '$\min_{s\in[t_s,T]} A(s)=0$');

lgd1 = legend(ax1, [hNo, hPartial, hFull, hSafety1], ...
    {'No major upgrade', 'Partial upgrade', 'Full upgrade', ...
     '$\min_{s\in[t_s,T]} A(s)=0$'}, ...
    'Interpreter', 'latex', 'Location', 'northwest');
lgd1.Box = 'off';

% ---------------- Panel (b): deviation region --------------------------
ax2 = nexttile;

contourf(ax2, RGrid, AGrid, double(deviationRegion), ...
    [-0.5, 0.5, 1.5], 'LineColor', 'none');

colormap(ax2, [0.96 0.96 0.96; ...
               0.38 0.38 0.38]);

hold(ax2, 'on');

% Library-based upgrade boundary
[~, hEconomic] = contour(ax2, RGrid, AGrid, ...
    double(upgradeSelected), [0.5 0.5], ...
    'k--', 'LineWidth', 1.8);

% Dynamic carbon-account safety boundary
[~, hSafety2] = contour(ax2, RGrid, AGrid, bestMinA, ...
    [0 0], 'k-', 'LineWidth', 2.0);

hold(ax2, 'off');

xlabel(ax2, 'Privacy risk, $R(t_s)$', 'Interpreter', 'latex');
ylabel(ax2, 'Carbon-account balance, $A(t_s^-)$', 'Interpreter', 'latex');
title(ax2, '(b) Upgrade selection and account deficit');
xlim(ax2, [RGrid(1), RGrid(end)]);
ylim(ax2, [AGrid(1), AGrid(end)]);
box(ax2, 'on');
grid(ax2, 'on');
set(ax2, 'Layer', 'top', 'LineWidth', 1.0);

hDeviation = patch(ax2, NaN, NaN, [0.38 0.38 0.38], ...
    'EdgeColor', 'k', ...
    'DisplayName', 'Selected upgrading with $A(s)<0$');
set(hEconomic, 'DisplayName', 'Library-based upgrade boundary');
set(hSafety2, 'DisplayName', '$\min_{s\in[t_s,T]} A(s)=0$');

lgd2 = legend(ax2, [hDeviation, hEconomic, hSafety2], ...
    {'Selected upgrading with $A(s)<0$', ...
     'Library-based upgrade boundary', '$\min_{s\in[t_s,T]} A(s)=0$'}, ...
    'Interpreter', 'latex', 'Location', 'northwest');
lgd2.Box = 'off';

%% 7. Numerical summary and optional export

fprintf('\n================ Figure 6 diagnostic =================\n');
fprintf('Decision slice: t = %.2f; remaining horizon = %.2f\n', tSlice, H);
fprintf('Fixed efficiency state on the slice: kappa = %.2f\n', kappaSlice);
fprintf('Number of R-A grid points: %d\n', numel(RR));
fprintf('Number of candidate policies: %d\n', nCandidates);
fprintf('All plotted states satisfy A(t_s) >= 0.\n');
fprintf('No-upgrade points: %d\n', sum(upgradeRegimePlot(:) == 0));
fprintf('Partial-upgrade points: %d\n', sum(upgradeRegimePlot(:) == 1));
fprintf('Full-upgrade points: %d\n', sum(upgradeRegimePlot(:) == 2));
fprintf('Dynamically unsafe points: %d\n', sum(~dynamicSafe(:)));
fprintf('Economic-safety deviation points: %d\n', sum(deviationRegion(:)));

% Diagnostic: the figure itself DOES NOT establish Proposition 5's triple
% separation unless selected upgrades also satisfy BAccount>0.
paybackPositive = upgradeSelected & (bestBAccount > 1e-10);
jumpDeficit = upgradeSelected & (bestUpgradeAplus < -1e-10);
fprintf('Selected upgrades with positive account payback: %d\n', ...
    sum(paybackPositive(:)));
fprintf('Unsafe selected upgrades with positive payback: %d\n', ...
    sum(sum(deviationRegion & paybackPositive)));
fprintf('Selected upgrades creating an immediate jump deficit: %d\n', ...
    sum(jumpDeficit(:)));
fprintf('Unsafe selected upgrades not attributable to immediate jump: %d\n', ...
    sum(sum(deviationRegion & ~jumpDeficit)));
fprintf('Elapsed time: %.2f seconds\n', toc);

% Set true when final files should be exported to the current folder.
saveOutput = false;

if saveOutput
    savefig(fig, 'Figure6_EconomicCarbonSafety.fig');
    exportgraphics(fig, 'Figure6_EconomicCarbonSafety.pdf', ...
        'ContentType', 'vector');
    exportgraphics(fig, 'Figure6_EconomicCarbonSafety.png', ...
        'Resolution', 600);
end
