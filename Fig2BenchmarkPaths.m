function Figure2_BenchmarkPaths_Final
%% ========================================================================
%  Fast standalone code for manuscript Figure 2
%  Benchmark dynamic paths under the dual policy
%
%  Speed-up:
%    - 21 maintenance decision knots instead of 121 decision variables;
%    - one coarse solve for each L=0,...,Lmax;
%    - three refinement starts only for the best L;
%    - normally 8 fmincon calls rather than about 50.
%
%  Final carbon-account model:
%    dA/dt = Q - chi*y/kappa - em*m - 0.5*em2*m^2
%    Phi_A = vA*A(T) - 0.5*phiA*max(-A(T),0)^2
%
%  Save this file exactly as Figure2_BenchmarkPaths_Final.m, then click Run.
%% ========================================================================

clc;
close all;
rng(20260910, 'twister');

p = baseline_parameters();
validate_parameters(p);

tSim = linspace(0, p.T, p.Nsim + 1)';
tM   = linspace(0, p.T, p.Nm)';

if exist('fmincon', 'file') ~= 2
    error('Optimization Toolbox is required because this code uses fmincon.');
end

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', p.solverDisplay, ...
    'MaxIterations', 600, ...
    'MaxFunctionEvaluations', 1.5e5, ...
    'StepTolerance', 1e-7, ...
    'OptimalityTolerance', 1e-5, ...
    'ConstraintTolerance', 1e-6);

%% ========================== Fast dual-policy solve ======================

fprintf('\nStage 1: coarse search over the number of upgrades ...\n');

coarseResults = cell(p.Lmax + 1, 1);
bestCoarseValue = -Inf;
bestL = 0;

for L = 0:p.Lmax
    coarseResults{L+1} = solve_fixed_L( ...
        L, p, tSim, tM, options, 1, []);

    fprintf('  L=%d, objective=%12.6f\n', L, coarseResults{L+1}.J);

    if coarseResults{L+1}.J > bestCoarseValue
        bestCoarseValue = coarseResults{L+1}.J;
        bestL = L;
    end
end

fprintf('\nStage 2: refine the best upgrade count L=%d ...\n', bestL);

res_dual = solve_fixed_L( ...
    bestL, p, tSim, tM, options, p.nRefineStarts, ...
    coarseResults{bestL+1}.z);

% The refinement contains the coarse solution as its first initial point,
% so it cannot be worse than the selected coarse result.
fprintf('Refined objective value = %.8f\n', res_dual.J);

%% ========================== Draw manuscript Figure 2 ===================

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 11);

upgradeTime  = res_dual.sigma(:);
upgradeOmega = res_dual.omega(:);
upgradeSize  = 80 + 140*upgradeOmega;

if isempty(upgradeTime)
    mAtUpgrade = [];
else
    mAtUpgrade = interp1(tSim, res_dual.mGrid, upgradeTime, ...
        'linear', 'extrap');
end

figure('Name', 'Figure 2: Benchmark optimal paths', ...
       'Color', 'w', 'Position', [100, 60, 860, 900]);

tl = tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

% ------------------------- (a) Maintenance effort ------------------------
ax1 = nexttile;
plot(ax1, tSim, res_dual.mGrid, 'k', 'LineWidth', 1.8);
hold(ax1, 'on');

for i = 1:numel(upgradeTime)
    xline(ax1, upgradeTime(i), ':', ...
        'Color', [0.70 0.70 0.70], ...
        'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
end

if ~isempty(upgradeTime)
    scatter(ax1, upgradeTime, mAtUpgrade, upgradeSize, 'o', ...
        'MarkerFaceColor', [0.25 0.25 0.25], ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.8);
end

hold(ax1, 'off');
ylabel(ax1, '$m^*(t)$', 'Interpreter', 'latex');
title(ax1, '(a) Optimal maintenance effort');
xlim(ax1, [0 p.T]);
ylim(ax1, [0 max(1.05*p.mbar, 1.05*max(res_dual.mGrid)+1e-3)]);
grid(ax1, 'on');
box(ax1, 'on');
set(ax1, 'LineWidth', 1.0);

% ---------------------- (b) Risk and efficiency --------------------------
ax2 = nexttile;

yyaxis(ax2, 'left');
hR = plot(ax2, res_dual.pathTime, res_dual.Rpath, ...
    'k', 'LineWidth', 1.8);
ylabel(ax2, '$R(t)$', 'Interpreter', 'latex');
ax2.YColor = [0 0 0];

yyaxis(ax2, 'right');
hK = plot(ax2, res_dual.pathTime, res_dual.Kpath, '--', ...
    'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.8);
ylabel(ax2, '$\kappa(t)$', 'Interpreter', 'latex');
ax2.YColor = [0.35 0.35 0.35];

hold(ax2, 'on');
for i = 1:numel(upgradeTime)
    xline(ax2, upgradeTime(i), ':', ...
        'Color', [0.70 0.70 0.70], ...
        'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
end
hold(ax2, 'off');

title(ax2, '(b) Privacy risk and efficiency dynamics');
xlim(ax2, [0 p.T]);
grid(ax2, 'on');
box(ax2, 'on');
set(ax2, 'LineWidth', 1.0);

lgd = legend(ax2, [hR, hK], {'$R(t)$', '$\kappa(t)$'}, ...
    'Interpreter', 'latex', ...
    'Location', 'best');
lgd.Box = 'off';
lgd.AutoUpdate = 'off';

% -------------------------- (c) Carbon account ---------------------------
ax3 = nexttile;
plot(ax3, res_dual.pathTime, res_dual.Apath, ...
    'k', 'LineWidth', 2.0);
hold(ax3, 'on');

yline(ax3, 0, '--', ...
    'Color', [0.45 0.45 0.45], ...
    'LineWidth', 1.2);

for i = 1:numel(upgradeTime)
    xline(ax3, upgradeTime(i), ':', ...
        'Color', [0.65 0.65 0.65], ...
        'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
end

hold(ax3, 'off');
ylabel(ax3, '$A(t)$', 'Interpreter', 'latex');
xlabel(ax3, 'Time');
title(ax3, '(c) Carbon account trajectory');
xlim(ax3, [0 p.T]);
grid(ax3, 'on');
box(ax3, 'on');
set(ax3, 'LineWidth', 1.0);

title(tl, 'Benchmark optimal paths under the dual policy', ...
    'FontWeight', 'normal');

% Optional export
% exportgraphics(gcf, 'Figure2_BenchmarkPaths.pdf', ...
%     'ContentType', 'vector');
% exportgraphics(gcf, 'Figure2_BenchmarkPaths.png', ...
%     'Resolution', 600);

%% ============================== Print results ===========================

fprintf('\n================ Figure 2 optimization result ================\n');
fprintf('Optimal number of active upgrades L = %d\n', res_dual.L);
fprintf('Objective value                    = %.8f\n', res_dual.J);
fprintf('Minimum carbon-account balance     = %.8f\n', res_dual.minA);
fprintf('Terminal carbon-account balance    = %.8f\n', res_dual.terminalA);
fprintf('Integrated maintenance effort      = %.8f\n', res_dual.maintenanceCum);
fprintf('First-order optimality             = %.3e\n', res_dual.firstOrderOpt);
fprintf('Constraint violation               = %.3e\n', ...
    res_dual.constraintViolation);

if isempty(res_dual.sigma)
    fprintf('Active upgrades: none\n');
else
    fprintf('Upgrade times:       ');
    fprintf('%7.3f ', res_dual.sigma);
    fprintf('\n');

    fprintf('Upgrade intensities: ');
    fprintf('%7.3f ', res_dual.omega);
    fprintf('\n');
end

fprintf(['\nThe baseline model permits temporary negative values of A(t); ' ...
         'A(t)>=0 is not imposed here.\n']);

end


%% ========================================================================
%                              Parameters
%% ========================================================================

function p = baseline_parameters()

p.T    = 10;
p.Nsim = 120;   % state-simulation intervals
p.Nm   = 21;    % maintenance decision knots; use 31 for a finer final check
p.r    = 0.05;

p.y     = 1.00;
p.price = 3.00;
p.pE    = 0.60;
p.chi   = 0.35;
p.Q     = 0.50;

p.R0       = 0.40;
p.kappa0   = 0.55;
p.A0       = 2.20;
p.kappaBar = 1.00;

p.beta   = 0.25;
p.deltaR = 0.05;
p.etaR   = 0.45;

p.etaK   = 0.35;
p.deltaK = 0.12;

p.cm  = 0.08;
p.am  = 0.45;
p.em  = 0.03;
p.em2 = 0.10;

p.aR = 1.20;
p.ay = 0.45;

% Upgrade economic cost: K0 + K1*omega + 0.5*K2*omega^2.
% Align these values with the manuscript parameter table.
p.K0 = 0.08;
p.K1 = 0.25;
p.K2 = 0.70;

% One-off upgrade carbon emissions.
p.xi1 = 0.14;
p.xi2 = 0.35;

p.gammaR = 0.65;
p.gammaK = 0.70;

% Terminal values.
p.vA   = 0.50;
p.phiA = 8.00;
p.vK   = 0.80;
p.phiR = 1.50;

p.mbar = 1.00;

% Endogenous upgrade settings.
p.Lmax     = 4;
p.sigmaMin = 0.50;
p.sigmaMax = 9.50;
p.minGap   = 0.80;
p.omegaMin = 0.05;

% Fast-search settings.
p.nRefineStarts = 3;
p.solverDisplay = 'off';  % use 'iter' only for diagnostics

end


function validate_parameters(p)

if p.Nsim < 2 || p.Nm < 2 || p.T <= 0
    error('Nsim and Nm must be at least 2, and T must be positive.');
end

if p.kappa0 <= 0 || p.kappa0 > p.kappaBar
    error('kappa0 must lie in (0,kappaBar].');
end

if p.omegaMin <= 0 || p.omegaMin > 1
    error('omegaMin must lie in (0,1].');
end

if p.sigmaMin < 0 || p.sigmaMax > p.T || ...
        p.sigmaMin >= p.sigmaMax
    error('Require 0<=sigmaMin<sigmaMax<=T.');
end

if p.Lmax >= 2
    if (p.Lmax-1)*p.minGap > p.sigmaMax-p.sigmaMin+1e-12
        error('Lmax and minGap are infeasible within the time bounds.');
    end
end

end


%% ========================================================================
%                              Optimization
%% ========================================================================

function res = solve_fixed_L(L, p, tSim, tM, options, nStarts, seedZ)

nM = p.Nm;
nvar = nM + 2*L;

idxM = 1:nM;
idxSigma = nM + (1:L);
idxOmega = nM + L + (1:L);

lb = zeros(nvar, 1);
ub = zeros(nvar, 1);

lb(idxM) = 0;
ub(idxM) = p.mbar;

if L > 0
    lb(idxSigma) = p.sigmaMin;
    ub(idxSigma) = p.sigmaMax;
    lb(idxOmega) = p.omegaMin;
    ub(idxOmega) = 1;
end

[Aineq, bineq] = spacing_constraints(nvar, idxSigma, p.minGap);
z0Set = initial_points(L, p, tM, nStarts, seedZ);

bestObj = Inf;
bestZ = [];
bestExitflag = -Inf;
bestOutput = struct();

objFun = @(z) objective_wrapper(z, L, p, tSim, tM);

for k = 1:size(z0Set, 2)
    z0 = z0Set(:, k);

    [zSol, objVal, exitflag, output] = fmincon( ...
        objFun, z0, Aineq, bineq, [], [], lb, ub, [], options);

    if isfinite(objVal) && objVal < bestObj
        bestObj = objVal;
        bestZ = zSol;
        bestExitflag = exitflag;
        bestOutput = output;
    end
end

if isempty(bestZ)
    error('No finite solution was obtained for L=%d.', L);
end

[mKnots, sigma, omega] = unpack_controls(bestZ, L, p);
res = simulate_policy(mKnots, sigma, omega, p, tSim, tM);

res.L = L;
res.z = bestZ;
res.exitflag = bestExitflag;

if isfield(bestOutput, 'firstorderopt')
    res.firstOrderOpt = bestOutput.firstorderopt;
else
    res.firstOrderOpt = NaN;
end

if isfield(bestOutput, 'constrviolation')
    res.constraintViolation = bestOutput.constrviolation;
else
    res.constraintViolation = max([0; Aineq*bestZ-bineq]);
end

end


function obj = objective_wrapper(z, L, p, tSim, tM)

[mKnots, sigma, omega] = unpack_controls(z, L, p);
res = simulate_policy(mKnots, sigma, omega, p, tSim, tM);

if isfinite(res.J)
    obj = -res.J;
else
    obj = 1e20;
end

end


function [mKnots, sigma, omega] = unpack_controls(z, L, p)

nM = p.Nm;
mKnots = z(1:nM);

if L == 0
    sigma = [];
    omega = [];
else
    sigma = z(nM + (1:L));
    omega = z(nM + L + (1:L));
end

end


function [Aineq, bineq] = spacing_constraints(nvar, idxSigma, minGap)

L = numel(idxSigma);

if L <= 1
    Aineq = zeros(0, nvar);
    bineq = zeros(0, 1);
    return;
end

Aineq = zeros(L-1, nvar);
bineq = -minGap*ones(L-1, 1);

for i = 1:L-1
    Aineq(i, idxSigma(i)) = 1;
    Aineq(i, idxSigma(i+1)) = -1;
end

end


function z0Set = initial_points(L, p, tM, nStarts, seedZ)

nM = p.Nm;
nvar = nM + 2*L;
z0Set = zeros(nvar, nStarts);

for k = 1:nStarts
    phase = 2*pi*(k-1)/max(nStarts,1);
    level = 0.35 + 0.35*(k-1)/max(nStarts-1,1);

    m0 = level + 0.05*sin(2*pi*tM/p.T + phase);
    m0 = min(max(m0, 0), p.mbar);

    z0 = zeros(nvar, 1);
    z0(1:nM) = m0;

    if L > 0
        if k == 1
            if L == 1
                sigma0 = 0.5*(p.sigmaMin+p.sigmaMax);
            else
                sigma0 = linspace(p.sigmaMin, p.sigmaMax, L)';
            end
        else
            sigma0 = random_feasible_times(L, p);
        end

        omegaLevel = 0.45 + 0.25*(k-1)/max(nStarts-1,1);
        omega0 = omegaLevel*ones(L, 1);
        omega0 = min(max(omega0, p.omegaMin), 1);

        z0(nM+(1:L)) = sigma0;
        z0(nM+L+(1:L)) = omega0;
    end

    z0Set(:, k) = z0;
end

if ~isempty(seedZ)
    if numel(seedZ) ~= nvar
        error('The refinement seed has an incorrect dimension.');
    end
    z0Set(:,1) = seedZ;
end

end


function sigma = random_feasible_times(L, p)

if L == 1
    sigma = p.sigmaMin + (p.sigmaMax-p.sigmaMin)*rand;
    return;
end

base = p.sigmaMin + (0:L-1)'*p.minGap;
remainingSpan = p.sigmaMax-base(end);
sigma = base + sort(rand(L,1))*remainingSpan;

end


%% ========================================================================
%                              State simulation
%% ========================================================================

function res = simulate_policy(mKnots, sigma, omega, p, tSim, tM)

mGrid = interp1(tM, mKnots, tSim, 'linear');
mGrid = min(max(mGrid, 0), p.mbar);

allTimes = sort([tSim(:); sigma(:)]);
tol = 1e-9;
timeList = allTimes([true; diff(allTimes)>tol]);

R = p.R0;
K = p.kappa0;
A = p.A0;

pathTime = 0;
Rpath = R;
Kpath = K;
Apath = A;

runningProfit = 0;
runningCarbon = 0;
riskCum = 0;
kappaCum = 0;
maintenanceCum = 0;
upgradeCostPV = 0;
upgradeCarbon = 0;

eventApplied = false(numel(sigma), 1);

for q = 1:numel(timeList)-1

    t0 = timeList(q);
    t1 = timeList(q+1);
    dt = t1-t0;

    if dt > 0
        tm = 0.5*(t0+t1);
        mMid = interp1(tM, mKnots, tm, 'linear', 'extrap');
        mMid = min(max(mMid,0),p.mbar);

        R0 = R;
        K0state = K;

        [R,K,A] = propagate_state(R,K,A,mMid,dt,p);

        Ravg = 0.5*(R0+R);
        Kavg = 0.5*(K0state+K);
        kEff = max(Kavg,1e-6);

        operatingCost = p.pE*p.y/kEff ...
                      + 0.5*p.ay*(p.y/kEff)^2;
        maintenanceCost = p.cm*mMid+0.5*p.am*mMid^2;
        privacyCost = 0.5*p.aR*Ravg^2;

        currentProfit = p.price*p.y ...
                      - operatingCost ...
                      - maintenanceCost ...
                      - privacyCost;

        runningProfit = runningProfit ...
                      + exp(-p.r*tm)*dt*currentProfit;

        currentCarbon = p.chi*p.y/kEff ...
                      + p.em*mMid ...
                      + 0.5*p.em2*mMid^2;

        runningCarbon = runningCarbon+dt*currentCarbon;
        riskCum = riskCum+dt*Ravg;
        kappaCum = kappaCum+dt*Kavg;
        maintenanceCum = maintenanceCum+dt*mMid;

        pathTime(end+1,1) = t1; %#ok<AGROW>
        Rpath(end+1,1) = R; %#ok<AGROW>
        Kpath(end+1,1) = K; %#ok<AGROW>
        Apath(end+1,1) = A; %#ok<AGROW>
    end

    idxEvents = find(abs(sigma-t1)<1e-7 & ~eventApplied);

    for jj = 1:numel(idxEvents)
        l = idxEvents(jj);
        w = omega(l);

        K = K+p.gammaK*w*(p.kappaBar-K);
        R = (1-p.gammaR*w)*R;

        eventCarbon = p.xi1*w+0.5*p.xi2*w^2;
        eventCost = p.K0+p.K1*w+0.5*p.K2*w^2;

        A = A-eventCarbon;
        upgradeCarbon = upgradeCarbon+eventCarbon;
        upgradeCostPV = upgradeCostPV ...
                      + exp(-p.r*sigma(l))*eventCost;

        pathTime(end+1,1) = t1; %#ok<AGROW>
        Rpath(end+1,1) = R; %#ok<AGROW>
        Kpath(end+1,1) = K; %#ok<AGROW>
        Apath(end+1,1) = A; %#ok<AGROW>

        eventApplied(l) = true;
    end
end

Adeficit = max(-A,0);

terminalValue = p.vA*A ...
              - 0.5*p.phiA*Adeficit^2 ...
              + p.vK*K ...
              - 0.5*p.phiR*R^2;

J = runningProfit ...
  - upgradeCostPV ...
  + exp(-p.r*p.T)*terminalValue;

res.J = J;
res.mKnots = mKnots;
res.mGrid = mGrid;
res.sigma = sigma;
res.omega = omega;

res.pathTime = pathTime;
res.Rpath = Rpath;
res.Kpath = Kpath;
res.Apath = Apath;

res.riskCum = riskCum;
res.kappaAvg = kappaCum/p.T;
res.carbonCum = runningCarbon+upgradeCarbon;
res.maintenanceCum = maintenanceCum;
res.minA = min(Apath);
res.terminalA = A;

end


function [Rnext,Knext,Anext] = propagate_state(R,K,A,m,dt,p)

x = [R;K;A];
f = @(xv) state_rhs(xv,m,p);

k1 = f(x);
k2 = f(x+0.5*dt*k1);
k3 = f(x+0.5*dt*k2);
k4 = f(x+dt*k3);

xNext = x+dt*(k1+2*k2+2*k3+k4)/6;

Rnext = max(xNext(1),0);
Knext = min(max(xNext(2),1e-6),p.kappaBar);
Anext = xNext(3);

end


function dx = state_rhs(x,m,p)

R = max(x(1),0);
K = min(max(x(2),1e-6),p.kappaBar);

dR = p.beta*p.y-(p.deltaR+p.etaR*m)*R;
dK = p.etaK*m*(p.kappaBar-K)-p.deltaK*K;

% Final baseline carbon account: no deltaA term.
dA = p.Q ...
   - p.chi*p.y/K ...
   - p.em*m ...
   - 0.5*p.em2*m^2;

dx = [dR;dK;dA];

end
