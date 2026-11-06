function Section6_2_HardCarbonConstraint_Simplified
% SECTION6_2_HARDCARBONCONSTRAINT_SIMPLIFIED
% Fast, publication-oriented simulation for Section 6.2.
%
% This version removes the full A0 sweep and retains two representative
% states:
%   1) a dynamically safe state, used to verify that the hard constraint
%      is nonbinding and the compliance cost is zero;
%   2) an economic-carbon-safety divergence state, used to identify the
%      strategy adjustment and compliance cost after A(t)>=0 is imposed.
%
% The economic model is not simplified. To reduce runtime, maintenance is
% parameterized by a small number of control nodes and linearly interpolated
% on the full state-simulation grid. The code still optimizes:
%   - the continuous maintenance path m(t);
%   - the endogenous upgrade count L=0,...,Lmax;
%   - upgrade times sigma_l;
%   - upgrade intensities omega_l.
%
% Main outputs
%   Figure_Section6_2_Simplified.pdf/png
%   Table_Section6_2_ScenarioComparison.csv
%   Table_Section6_2_BindingChanges.csv
%
% Save this file with the exact name
%   Section6_2_HardCarbonConstraint_Simplified.m
% and run
%   Section6_2_HardCarbonConstraint_Simplified
%
% Required toolbox: MATLAB Optimization Toolbox (fmincon).

clearvars; clc; close all;
rng(20260918, 'twister');

p = baseline_parameters();
validate_parameters(p);

if exist('fmincon', 'file') ~= 2
    error('This program requires fmincon from MATLAB Optimization Toolbox.');
end

tGrid = linspace(0,p.T,p.N+1)';

options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display',p.solverDisplay, ...
    'MaxIterations',500, ...
    'MaxFunctionEvaluations',8e4, ...
    'StepTolerance',1e-8, ...
    'OptimalityTolerance',1e-6, ...
    'ConstraintTolerance',1e-7);

fprintf('\n===============================================================\n');
fprintf('Section 6.2: simplified hard-carbon-account simulation\n');
fprintf('===============================================================\n');
fprintf(['State nodes=%d, maintenance-control nodes=%d, ' ...
    'multistarts=%d, L=0,...,%d\n'], ...
    p.N+1,p.Ncontrol,p.nStarts,p.Lmax);

%% 1. Dynamically safe state

pSafe = p;
pSafe.A0 = p.A0Safe;

if p.runSafeCheck
    fprintf('\n--- Safe state: A0 = %.4f ---\n',pSafe.A0);
    resSafeSoft = solve_dual_model(pSafe,tGrid,options,false,[]);

    if resSafeSoft.minA < -p.hardAcceptTol
        error(['The chosen safe-state value A0Safe=%.4f is not dynamically ' ...
            'safe under the soft optimum (minA=%.6f). Increase p.A0Safe.'], ...
            pSafe.A0,resSafeSoft.minA);
    end

    % Feasible-set nesting implies that a dynamically safe unrestricted
    % optimum is also the hard optimum, so it is not solved a second time.
    resSafeHard = resSafeSoft;
    fprintf(['The soft optimum is dynamically safe; therefore the hard ' ...
        'constraint is nonbinding and Ccomp=0.\n']);
else
    fprintf(['\nSafe-state optimization skipped. The previous full-grid run ' ...
        'already verified min A>0 at A0=%.2f.\n'],pSafe.A0);
end

%% 2. Economic-carbon-safety divergence state

pBind = p;
pBind.A0 = p.A0Binding;

fprintf('\n--- Binding state: A0 = %.4f ---\n',pBind.A0);
resBindSoft = solve_dual_model(pBind,tGrid,options,false,[]);

if resBindSoft.minA >= -p.bindingTol
    error(['The chosen binding-state value A0Binding=%.4f does not create ' ...
        'a soft-model deficit (minA=%.6f). Reduce p.A0Binding or select the ' ...
        'economic-safety divergence point identified in Figure 6.'], ...
        pBind.A0,resBindSoft.minA);
end

% Warm-start the hard problem from the soft solution. Other multistarts are
% retained, so the hard strategy is fully reoptimized rather than truncated.
resBindHard = solve_dual_model(pBind,tGrid,options,true,resBindSoft);

if ~resBindHard.feasible
    error(['No feasible hard solution was obtained at A0Binding=%.4f. ' ...
        'Increase A0Binding slightly or raise p.nStarts.'],pBind.A0);
end

% Conditional nesting diagnostic. A soft optimum cannot theoretically be
% lower than the hard optimum. Re-polish only when numerical local optima
% violate this ordering; no repeated solve is used in the normal case.
if resBindHard.J > resBindSoft.J+p.valueTol
    warning(['The first run returned V^H>V^B. Re-solving the soft model ' ...
        'once from the hard solution to correct a local-optimum artifact.']);
    resSoftCheck = solve_dual_model(pBind,tGrid,options,false,resBindHard);
    if resSoftCheck.J > resBindSoft.J
        resBindSoft = resSoftCheck;
    end
end

CBind = resBindSoft.J-resBindHard.J;
PBind = 100*CBind/max(abs(resBindSoft.J),1e-10);

if CBind < -p.valueTol
    warning(['V^H still exceeds V^B by %.3e. Increase p.nStarts before ' ...
        'using the result in the manuscript.'],-CBind);
end

%% 3. Strategy and performance indicators

MBindSoft = trapz(tGrid,resBindSoft.mGrid);
MBindHard = trapz(tGrid,resBindHard.mGrid);
OmegaBindSoft = sum(resBindSoft.omega);
OmegaBindHard = sum(resBindHard.omega);
SigmaBindSoft = weighted_upgrade_time(resBindSoft.sigma,resBindSoft.omega);
SigmaBindHard = weighted_upgrade_time(resBindHard.sigma,resBindHard.omega);

DeltaM = MBindHard-MBindSoft;
DeltaOmega = OmegaBindHard-OmegaBindSoft;
DeltaSigma = SigmaBindHard-SigmaBindSoft;
DeltaRisk = resBindHard.riskCum-resBindSoft.riskCum;
DeltaKappa = resBindHard.kappaAvg-resBindSoft.kappaAvg;

PctM = relative_change(MBindHard,MBindSoft);
PctOmega = relative_change(OmegaBindHard,OmegaBindSoft);
PctSigma = relative_change(SigmaBindHard,SigmaBindSoft);
PctRisk = relative_change(resBindHard.riskCum,resBindSoft.riskCum);
PctKappa = relative_change(resBindHard.kappaAvg,resBindSoft.kappaAvg);

%% 4. Export result tables

if p.runSafeCheck
    MSafeSoft = trapz(tGrid,resSafeSoft.mGrid);
    MSafeHard = trapz(tGrid,resSafeHard.mGrid);
    OmegaSafeSoft = sum(resSafeSoft.omega);
    OmegaSafeHard = sum(resSafeHard.omega);
    SigmaSafeSoft = weighted_upgrade_time(resSafeSoft.sigma,resSafeSoft.omega);
    SigmaSafeHard = weighted_upgrade_time(resSafeHard.sigma,resSafeHard.omega);

    Scenario = {'Safe';'Safe';'Binding';'Binding'};
    Constraint = {'Soft';'Hard';'Soft';'Hard'};
    A0 = [pSafe.A0;pSafe.A0;pBind.A0;pBind.A0];
    ObjectiveValue = [resSafeSoft.J;resSafeHard.J;resBindSoft.J;resBindHard.J];
    MinimumAccount = [resSafeSoft.minA;resSafeHard.minA; ...
                      resBindSoft.minA;resBindHard.minA];
    TerminalAccount = [resSafeSoft.terminalA;resSafeHard.terminalA; ...
                       resBindSoft.terminalA;resBindHard.terminalA];
    UpgradeCount = [resSafeSoft.L;resSafeHard.L;resBindSoft.L;resBindHard.L];
    MaintenanceTotal = [MSafeSoft;MSafeHard;MBindSoft;MBindHard];
    UpgradeIntensityTotal = [OmegaSafeSoft;OmegaSafeHard; ...
                             OmegaBindSoft;OmegaBindHard];
    WeightedUpgradeTime = [SigmaSafeSoft;SigmaSafeHard; ...
                           SigmaBindSoft;SigmaBindHard];
    CumulativePrivacyRisk = [resSafeSoft.riskCum;resSafeHard.riskCum; ...
                             resBindSoft.riskCum;resBindHard.riskCum];
    AverageEfficiency = [resSafeSoft.kappaAvg;resSafeHard.kappaAvg; ...
                         resBindSoft.kappaAvg;resBindHard.kappaAvg];
    ComplianceCost = [0;0;CBind;CBind];
    RelativeComplianceCost = [0;0;PBind;PBind];
else
    Scenario = {'Binding';'Binding'};
    Constraint = {'Soft';'Hard'};
    A0 = [pBind.A0;pBind.A0];
    ObjectiveValue = [resBindSoft.J;resBindHard.J];
    MinimumAccount = [resBindSoft.minA;resBindHard.minA];
    TerminalAccount = [resBindSoft.terminalA;resBindHard.terminalA];
    UpgradeCount = [resBindSoft.L;resBindHard.L];
    MaintenanceTotal = [MBindSoft;MBindHard];
    UpgradeIntensityTotal = [OmegaBindSoft;OmegaBindHard];
    WeightedUpgradeTime = [SigmaBindSoft;SigmaBindHard];
    CumulativePrivacyRisk = [resBindSoft.riskCum;resBindHard.riskCum];
    AverageEfficiency = [resBindSoft.kappaAvg;resBindHard.kappaAvg];
    ComplianceCost = [CBind;CBind];
    RelativeComplianceCost = [PBind;PBind];
end

scenarioTable = table(Scenario,Constraint,A0,ObjectiveValue,MinimumAccount, ...
    TerminalAccount,UpgradeCount,MaintenanceTotal,UpgradeIntensityTotal, ...
    WeightedUpgradeTime,CumulativePrivacyRisk,AverageEfficiency, ...
    ComplianceCost,RelativeComplianceCost);

Indicator = {'Compliance cost';'Maintenance';'Upgrade intensity'; ...
             'Weighted upgrade time';'Privacy risk';'Average efficiency'};
AbsoluteChange = [CBind;DeltaM;DeltaOmega;DeltaSigma;DeltaRisk;DeltaKappa];
RelativeChangePercent = [PBind;PctM;PctOmega;PctSigma;PctRisk;PctKappa];

changeTable = table(Indicator,AbsoluteChange,RelativeChangePercent);

writetable(scenarioTable,'Table_Section6_2_ScenarioComparison.csv');
writetable(changeTable,'Table_Section6_2_BindingChanges.csv');

%% 5. Three-panel publication figure

set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',11);

softColor = [0.05 0.05 0.05];
hardColor = [0.50 0.50 0.50];

fig = figure('Name','Section 6.2: hard carbon-account constraint', ...
    'Color','w','Position',[60,180,1520,450]);
tl = tiledlayout(fig,1,3,'Padding','compact','TileSpacing','compact');

% (a) Carbon-account paths in the binding state.
ax1 = nexttile(tl);
hASoft = plot(ax1,resBindSoft.pathTime,resBindSoft.Apath,'-', ...
    'Color',softColor,'LineWidth',1.9);
hold(ax1,'on');
hAHard = plot(ax1,resBindHard.pathTime,resBindHard.Apath,'--', ...
    'Color',hardColor,'LineWidth',2.1);
yline(ax1,0,':','Color',[0.25 0.25 0.25],'LineWidth',1.2, ...
    'HandleVisibility','off');
hold(ax1,'off');
xlim(ax1,[0 p.T]);
xlabel(ax1,'Time');
ylabel(ax1,'$A(t)$','Interpreter','latex');
title(ax1,'(a) Dynamic carbon-account paths');
legend(ax1,[hASoft,hAHard],{'Soft constraint','Hard constraint'}, ...
    'Location','best','Box','off');
grid(ax1,'on'); box(ax1,'on'); set(ax1,'LineWidth',1.0);

% (b) Maintenance paths and upgrade impulses in the binding state.
ax2 = nexttile(tl);
yyaxis(ax2,'left');
hmSoft = plot(ax2,tGrid,resBindSoft.mGrid,'-', ...
    'Color',softColor,'LineWidth',1.8);
hold(ax2,'on');
hmHard = plot(ax2,tGrid,resBindHard.mGrid,'--', ...
    'Color',hardColor,'LineWidth',2.0);
ylabel(ax2,'$m^*(t)$','Interpreter','latex');
ylim(ax2,[0,max(p.mbar,1.05*max([resBindSoft.mGrid;resBindHard.mGrid]))]);
ax2.YColor = [0 0 0];

yyaxis(ax2,'right');
if isempty(resBindSoft.sigma)
    huSoft = plot(ax2,nan,nan,'o','Color',softColor, ...
        'MarkerFaceColor','w','LineStyle','none');
else
    huSoft = stem(ax2,resBindSoft.sigma,resBindSoft.omega,'o', ...
        'Color',softColor,'LineWidth',1.1, ...
        'MarkerFaceColor','w','MarkerSize',5);
end
if isempty(resBindHard.sigma)
    huHard = plot(ax2,nan,nan,'s','Color',hardColor, ...
        'MarkerFaceColor',hardColor,'LineStyle','none');
else
    huHard = stem(ax2,resBindHard.sigma,resBindHard.omega,'s', ...
        'Color',hardColor,'LineWidth',1.1, ...
        'MarkerFaceColor',hardColor,'MarkerSize',5);
end
ylabel(ax2,'Upgrade intensity, $\omega_l$','Interpreter','latex');
ylim(ax2,[0 1.05]);
ax2.YColor = [0.35 0.35 0.35];
hold(ax2,'off');
xlim(ax2,[0 p.T]);
xlabel(ax2,'Time');
title(ax2,'(b) Governance-strategy adjustment');
legend(ax2,[hmSoft,hmHard,huSoft,huHard], ...
    {'Maintenance: soft','Maintenance: hard', ...
     'Upgrade: soft','Upgrade: hard'}, ...
    'Location','best','Box','off','FontSize',9);
grid(ax2,'on'); box(ax2,'on'); set(ax2,'LineWidth',1.0);

% (c) Relative changes caused by the hard constraint.
ax3 = nexttile(tl);
barValues = [PBind,PctM,PctOmega,PctSigma,PctRisk,PctKappa];
b = bar(ax3,barValues,'FaceColor','w','EdgeColor','k','LineWidth',1.2);
hold(ax3,'on');
yline(ax3,0,'-','Color',[0.30 0.30 0.30],'LineWidth',1.0, ...
    'HandleVisibility','off');
hold(ax3,'off');

set(ax3,'XTick',1:6, ...
    'XTickLabel',{'Compliance cost','Maintenance','Upgrade intensity', ...
                  'Upgrade timing','Privacy risk','Efficiency'}, ...
    'XTickLabelRotation',30,'LineWidth',1.0);
ylabel(ax3,'Relative change (\%)','Interpreter','latex');
title(ax3,'(c) Compliance cost and performance effects');
grid(ax3,'on'); box(ax3,'on');

drawnow;
yl = ylim(ax3);
span = max(yl(2)-yl(1),1);
offset = 0.025*span;
for j = 1:numel(barValues)
    if isnan(barValues(j))
        text(ax3,j,0,'n.a.','HorizontalAlignment','center', ...
            'VerticalAlignment','bottom','FontSize',8);
    elseif barValues(j)>=0
        text(ax3,j,barValues(j)+offset,sprintf('%.2f',barValues(j)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',8);
    else
        text(ax3,j,barValues(j)-offset,sprintf('%.2f',barValues(j)), ...
            'HorizontalAlignment','center','VerticalAlignment','top', ...
            'FontSize',8);
    end
end

% title(tl,sprintf(['Hard dynamic carbon-account constraint: ' ...
%     '$A_0=%.2f$, $C_{\\mathrm{comp}}=%.4f$'],pBind.A0,CBind), ...
%     'Interpreter','latex','FontWeight','normal');

exportgraphics(fig,'Figure_Section6_2_Simplified.pdf', ...
    'ContentType','vector');
exportgraphics(fig,'Figure_Section6_2_Simplified.png', ...
    'Resolution',600);

%% 6. Console report

if p.runSafeCheck
    fprintf('\n===============================================================\n');
    fprintf('Safe state: A0 = %.4f\n',pSafe.A0);
    fprintf('===============================================================\n');
    print_result('Soft/Hard',resSafeSoft,MSafeSoft, ...
        OmegaSafeSoft,SigmaSafeSoft);
    fprintf('Constraint status: nonbinding; Ccomp=0; Pcomp=0%%.\n');
end

fprintf('\n===============================================================\n');
fprintf('Binding state: A0 = %.4f\n',pBind.A0);
fprintf('===============================================================\n');
print_result('Soft constraint',resBindSoft,MBindSoft, ...
    OmegaBindSoft,SigmaBindSoft);
print_result('Hard constraint',resBindHard,MBindHard, ...
    OmegaBindHard,SigmaBindHard);

fprintf('Absolute compliance cost Ccomp       = %12.8f\n',CBind);
fprintf('Relative compliance cost Pcomp       = %12.6f %%\n',PBind);
fprintf('Change in maintenance Delta M        = %12.8f\n',DeltaM);
fprintf('Change in upgrade intensity Delta Om = %12.8f\n',DeltaOmega);
fprintf('Change in weighted time Delta sigma  = %12.8f\n',DeltaSigma);
fprintf('Change in cumulative risk            = %12.8f\n',DeltaRisk);
fprintf('Change in average efficiency         = %12.8f\n',DeltaKappa);

fprintf('\nFiles written:\n');
fprintf('  Figure_Section6_2_Simplified.pdf\n');
fprintf('  Figure_Section6_2_Simplified.png\n');
fprintf('  Table_Section6_2_ScenarioComparison.csv\n');
fprintf('  Table_Section6_2_BindingChanges.csv\n');

end


%% ======================================================================
% Parameters: identical economic calibration to the final baseline model
%% ======================================================================

function p = baseline_parameters()

p.T  = 10;
p.N  = 120;
p.dt = p.T/p.N;
p.Ncontrol = 16;
p.r  = 0.05;

p.y     = 1.00;
p.price = 3.00;
p.pE    = 0.60;
p.chi   = 0.35;
p.Q     = 0.50;

p.R0       = 0.40;
p.kappa0   = 0.55;
p.kappaBar = 1.00;

% Representative initial carbon-account balances.
% Replace A0Binding by the exact divergence-state coordinate used in
% Figure 6 if that figure adopts a different value.
p.A0Safe = 2.20;
p.A0Binding = 0.80;

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

p.K0 = 0.08;
p.K1 = 0.25;
p.K2 = 0.70;

p.xi1 = 0.14;
p.xi2 = 0.35;

p.gammaR = 0.65;
p.gammaK = 0.70;

p.vA   = 0.50;
p.phiA = 8.00;
p.vK   = 0.80;
p.phiR = 1.50;

p.mbar = 1.00;

p.Lmax     = 4;
p.sigmaMin = 0.50;
p.sigmaMax = 9.50;
p.minGap   = 0.80;
p.omegaMin = 0.05;

% The earlier full-grid run already established that A0Safe=2.20 is
% nonbinding. Keep this false for the fast Section 6.2 run. Set it to true
% only if the safe-state result must be reproduced from scratch.
p.runSafeCheck = false;

% Three multistarts are used with 16 maintenance-control nodes. State paths
% are still integrated on all N+1=121 simulation nodes.
p.nStarts = 3;
p.solverDisplay = 'off';

p.bindingTol = 1e-4;
p.hardAcceptTol = 5e-6;
p.valueTol = 1e-6;

end


function validate_parameters(p)

if p.N<2 || p.T<=0
    error('N must be at least 2 and T must be positive.');
end
if p.Ncontrol<2 || p.Ncontrol>p.N+1
    error('Ncontrol must lie between 2 and N+1.');
end
if p.kappa0<=0 || p.kappa0>p.kappaBar
    error('kappa0 must lie in (0,kappaBar].');
end
if p.A0Safe<0 || p.A0Binding<0
    error('The two initial carbon-account balances must be nonnegative.');
end
if p.omegaMin<=0 || p.omegaMin>1
    error('omegaMin must lie in (0,1].');
end
if p.sigmaMin<0 || p.sigmaMax>p.T || p.sigmaMin>=p.sigmaMax
    error('Upgrade times must satisfy 0<=sigmaMin<sigmaMax<=T.');
end
if (p.Lmax-1)*p.minGap>p.sigmaMax-p.sigmaMin+1e-12
    error('Lmax and minGap are infeasible within the time bounds.');
end

end


%% ======================================================================
% Strategy-level optimization
%% ======================================================================

function bestRes = solve_dual_model(p,tGrid,options,hardFlag,warmRes)

bestValue = -Inf;
bestRes = empty_result();

for L = 0:p.Lmax
    modelName = 'soft';
    if hardFlag
        modelName = 'hard';
    end
    fprintf('  solving %-4s L=%d ...\n',modelName,L);
    blockTimer = tic;
    currentRes = solve_fixed_count(L,p,tGrid,options,hardFlag,warmRes);

    if currentRes.feasible
        fprintf('  %-4s L=%d, V=%11.6f, minA=%9.5f, time=%7.1f s\n', ...
            modelName,L,currentRes.J,currentRes.minA,toc(blockTimer));

        if currentRes.J>bestValue
            bestValue = currentRes.J;
            bestRes = currentRes;
        end
    else
        fprintf('  %-4s L=%d, no feasible solution, time=%7.1f s\n', ...
            modelName,L,toc(blockTimer));
    end
end

if ~bestRes.feasible && ~hardFlag
    error('No finite soft-model solution was obtained.');
end

end


function res = solve_fixed_count(L,p,tGrid,options,hardFlag,warmRes)

nM = p.Ncontrol;
nvar = nM+2*L;
idxM = 1:nM;
idxSigma = nM+(1:L);
idxOmega = nM+L+(1:L);

lb = zeros(nvar,1);
ub = zeros(nvar,1);
lb(idxM) = 0;
ub(idxM) = p.mbar;

if L>0
    lb(idxSigma) = p.sigmaMin;
    ub(idxSigma) = p.sigmaMax;
    lb(idxOmega) = p.omegaMin;
    ub(idxOmega) = 1;
end

[Aineq,bineq] = spacing_constraints(nvar,idxSigma,p.minGap);
z0Set = build_initial_points(L,p,tGrid,nvar,idxM,idxSigma,idxOmega);

if ~isempty(warmRes) && isstruct(warmRes) && warmRes.feasible && warmRes.L==L
    zWarm = pack_controls(warmRes.mCtrl,warmRes.sigma,warmRes.omega);
    if numel(zWarm)==nvar
        z0Set = [zWarm(:),z0Set]; %#ok<AGROW>
    end
end

objFun = @(z) objective_wrapper(z,L,p,tGrid);
if hardFlag
    nonlcon = @(z) hard_path_constraint(z,L,p,tGrid);
else
    nonlcon = [];
end

bestObj = Inf;
bestZ = [];
bestExitflag = -Inf;
bestOutput = struct();
bestViolation = Inf;

for k = 1:size(z0Set,2)
    z0 = z0Set(:,k);
    startTimer = tic;
    fprintf('      multistart %d/%d ... ',k,size(z0Set,2));

    [zSol,objVal,exitflag,output] = fmincon(objFun,z0, ...
        Aineq,bineq,[],[],lb,ub,nonlcon,options);
    fprintf('%.1f s\n',toc(startTimer));

    if ~isfinite(objVal)
        continue;
    end

    [mCtrl,sigma,omega] = unpack_controls(zSol,L,p);
    trial = simulate_policy(mCtrl,sigma,omega,p,tGrid);
    pathViolation = max(0,-trial.minA);

    if ~hardFlag || pathViolation<=p.hardAcceptTol
        if objVal<bestObj
            bestObj = objVal;
            bestZ = zSol;
            bestExitflag = exitflag;
            bestOutput = output;
            bestViolation = pathViolation;
        end
    end
end

if isempty(bestZ)
    res = empty_result();
    res.L = L;
    return;
end

[mCtrl,sigma,omega] = unpack_controls(bestZ,L,p);
res = simulate_policy(mCtrl,sigma,omega,p,tGrid);
res.L = L;
res.exitflag = bestExitflag;
res.pathViolation = bestViolation;
res.feasible = isfinite(res.J) && (~hardFlag || bestViolation<=p.hardAcceptTol);

if isfield(bestOutput,'firstorderopt')
    res.firstOrderOpt = bestOutput.firstorderopt;
else
    res.firstOrderOpt = NaN;
end
if isfield(bestOutput,'constrviolation')
    res.solverConstraintViolation = bestOutput.constrviolation;
else
    res.solverConstraintViolation = NaN;
end

end


function obj = objective_wrapper(z,L,p,tGrid)

[mCtrl,sigma,omega] = unpack_controls(z,L,p);
res = simulate_policy(mCtrl,sigma,omega,p,tGrid);

if isfinite(res.J)
    obj = -res.J;
else
    obj = 1e20;
end

end


function [c,ceq] = hard_path_constraint(z,L,p,tGrid)

[mCtrl,sigma,omega] = unpack_controls(z,L,p);
res = simulate_policy(mCtrl,sigma,omega,p,tGrid);

% Fixed-dimension equivalent of A(t)>=0 over all simulated path nodes,
% including the balances immediately before and after every upgrade.
c = -res.minA;
ceq = [];

end


function [Aineq,bineq] = spacing_constraints(nvar,idxSigma,minGap)

L = numel(idxSigma);
if L<=1
    Aineq = zeros(0,nvar);
    bineq = zeros(0,1);
    return;
end

Aineq = zeros(L-1,nvar);
bineq = -minGap*ones(L-1,1);
for i = 1:L-1
    Aineq(i,idxSigma(i)) = 1;
    Aineq(i,idxSigma(i+1)) = -1;
end

end


function z0Set = build_initial_points(L,p,tGrid,nvar,idxM,idxSigma,idxOmega)

nStarts = p.nStarts;
z0Set = zeros(nvar,nStarts);
maintenanceLevels = linspace(0.15,0.85,nStarts);
omegaLevels = linspace(max(p.omegaMin,0.15),0.85,nStarts);
tControl = linspace(0,p.T,p.Ncontrol)';

for k = 1:nStarts
    z0 = zeros(nvar,1);
    phase = 2*pi*(k-1)/max(nStarts,1);
    m0 = maintenanceLevels(k)+0.05*sin(2*pi*tControl/p.T+phase);
    z0(idxM) = min(max(m0,0),p.mbar);

    if L>0
        if k==1
            if L==1
                sigma0 = 0.5*(p.sigmaMin+p.sigmaMax);
            else
                sigma0 = linspace(p.sigmaMin,p.sigmaMax,L)';
            end
        else
            sigma0 = random_feasible_times(L,p);
        end

        omega0 = omegaLevels(k)*ones(L,1);
        if k>2
            omega0 = omega0+0.08*randn(L,1);
        end
        omega0 = min(max(omega0,p.omegaMin),1);

        z0(idxSigma) = sigma0;
        z0(idxOmega) = omega0;
    end

    z0Set(:,k) = z0;
end

end


function sigma = random_feasible_times(L,p)

if L==0
    sigma = [];
elseif L==1
    sigma = p.sigmaMin+(p.sigmaMax-p.sigmaMin)*rand;
else
    base = p.sigmaMin+(0:L-1)'*p.minGap;
    remainingSpan = p.sigmaMax-base(end);
    sigma = base+sort(rand(L,1))*remainingSpan;
end

end


function [mCtrl,sigma,omega] = unpack_controls(z,L,p)

nM = p.Ncontrol;
mCtrl = z(1:nM);

if L>0
    sigma = z(nM+(1:L));
    omega = z(nM+L+(1:L));
else
    sigma = [];
    omega = [];
end

end


function z = pack_controls(mCtrl,sigma,omega)

z = [mCtrl(:);sigma(:);omega(:)];

end


%% ======================================================================
% State simulation: identical to the final baseline model
%% ======================================================================

function res = simulate_policy(mCtrl,sigma,omega,p,tGrid)

% Piecewise-linear maintenance parameterization. Only Ncontrol values are
% optimized, while all state dynamics and the hard carbon-account test are
% evaluated on the original N+1 simulation grid and at upgrade events.
tControl = linspace(0,p.T,p.Ncontrol)';
mGrid = interp1(tControl,mCtrl(:),tGrid,'linear');
mGrid = min(max(mGrid,0),p.mbar);

eventTimes = sigma(:);
allTimes = sort([tGrid(:);eventTimes]);
tol = 1e-9;
timeList = allTimes([true;diff(allTimes)>tol]);

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
upgradeCostPV = 0;
upgradeCarbon = 0;
eventApplied = false(numel(sigma),1);

for q = 1:numel(timeList)-1
    t0 = timeList(q);
    t1 = timeList(q+1);
    dt = t1-t0;

    if dt>0
        tm = 0.5*(t0+t1);
        mMid = interp1(tGrid,mGrid,tm,'linear','extrap');
        mMid = min(max(mMid,0),p.mbar);

        R0state = R;
        K0state = K;
        [R,K,A] = propagate_state(R,K,A,mMid,dt,p);

        kEff = max(0.5*(K0state+K),1e-6);
        Ravg = 0.5*(R0state+R);
        Kavg = 0.5*(K0state+K);

        operatingCost = p.pE*p.y/kEff+0.5*p.ay*(p.y/kEff)^2;
        maintenanceCost = p.cm*mMid+0.5*p.am*mMid^2;
        privacyCost = 0.5*p.aR*Ravg^2;
        currentProfit = p.price*p.y-operatingCost-maintenanceCost-privacyCost;

        runningProfit = runningProfit+exp(-p.r*tm)*dt*currentProfit;

        currentCarbon = p.chi*p.y/kEff+p.em*mMid+0.5*p.em2*mMid^2;
        runningCarbon = runningCarbon+dt*currentCarbon;
        riskCum = riskCum+dt*Ravg;
        kappaCum = kappaCum+dt*Kavg;

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
        A = A-eventCarbon;

        eventCost = p.K0+p.K1*w+0.5*p.K2*w^2;
        upgradeCostPV = upgradeCostPV+exp(-p.r*sigma(l))*eventCost;
        upgradeCarbon = upgradeCarbon+eventCarbon;

        pathTime(end+1,1) = t1; %#ok<AGROW>
        Rpath(end+1,1) = R; %#ok<AGROW>
        Kpath(end+1,1) = K; %#ok<AGROW>
        Apath(end+1,1) = A; %#ok<AGROW>
        eventApplied(l) = true;
    end
end

Adeficit = max(-A,0);
terminalValue = p.vA*A-0.5*p.phiA*Adeficit^2 ...
    +p.vK*K-0.5*p.phiR*R^2;
J = runningProfit-upgradeCostPV+exp(-p.r*p.T)*terminalValue;

res.mCtrl = mCtrl(:);
res.mGrid = mGrid;
res.sigma = sigma;
res.omega = omega;
res.pathTime = pathTime;
res.Rpath = Rpath;
res.Kpath = Kpath;
res.Apath = Apath;
res.J = J;
res.runningProfit = runningProfit;
res.upgradeCostPV = upgradeCostPV;
res.terminalValue = terminalValue;
res.riskCum = riskCum;
res.carbonCum = runningCarbon+upgradeCarbon;
res.runningCarbon = runningCarbon;
res.upgradeCarbon = upgradeCarbon;
res.kappaAvg = kappaCum/p.T;
res.minA = min(Apath);
res.terminalA = A;
res.terminalR = R;
res.terminalKappa = K;
res.feasible = isfinite(J);
res.L = numel(sigma);
res.exitflag = NaN;
res.firstOrderOpt = NaN;
res.solverConstraintViolation = NaN;
res.pathViolation = max(0,-res.minA);

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
dA = p.Q-p.chi*p.y/K-p.em*m-0.5*p.em2*m^2;

dx = [dR;dK;dA];

end


%% ======================================================================
% Reporting helpers
%% ======================================================================

function value = weighted_upgrade_time(sigma,omega)

if isempty(omega) || sum(omega)<=1e-12
    value = NaN;
else
    value = sum(sigma(:).*omega(:))/sum(omega);
end

end


function pct = relative_change(hardValue,softValue)

if ~isfinite(hardValue) || ~isfinite(softValue) || abs(softValue)<=1e-10
    pct = NaN;
else
    pct = 100*(hardValue-softValue)/abs(softValue);
end

end


function res = empty_result()

res = struct('feasible',false,'J',-Inf,'L',0,'minA',NaN, ...
    'terminalA',NaN,'mCtrl',[],'mGrid',[],'sigma',[],'omega',[], ...
    'pathTime',[],'Apath',[],'Rpath',[],'Kpath',[], ...
    'riskCum',NaN,'kappaAvg',NaN,'carbonCum',NaN, ...
    'pathViolation',Inf,'exitflag',NaN,'firstOrderOpt',NaN, ...
    'solverConstraintViolation',NaN);

end


function print_result(name,res,M,Omega,sigmaBar)

fprintf('\n%-18s V=%12.7f, minA=%10.6f, terminalA=%10.6f\n', ...
    name,res.J,res.minA,res.terminalA);
fprintf('%-18s L=%d, M=%10.6f, Omega=%10.6f, weighted sigma=%10.6f\n', ...
    name,res.L,M,Omega,sigmaBar);
fprintf('%-18s risk=%10.6f, average kappa=%10.6f\n', ...
    name,res.riskCum,res.kappaAvg);

if isempty(res.sigma)
    fprintf('%-18s upgrades: none\n',name);
else
    fprintf('%-18s sigma = ',name);
    fprintf('%7.3f ',res.sigma);
    fprintf('\n');
    fprintf('%-18s omega = ',name);
    fprintf('%7.3f ',res.omega);
    fprintf('\n');
end

fprintf('%-18s path violation=%9.3e, exitflag=%g\n', ...
    name,res.pathViolation,res.exitflag);

end
