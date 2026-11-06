function Section6_5_TraceDrivenCase
% SECTION6_5_TRACEDRIVENCASE
% Trace-driven representative case for the low-carbon computing model.
%
% The program keeps the final economic model unchanged and replaces the
% constant workload only through a time-aware workload interface. It solves:
%   1) constant workload, soft carbon-account model;
%   2) trace workload, soft carbon-account model;
%   3) trace workload, hard path constraint A(t)>=0.
%
% Set runMode='quick' for diagnostics. After all diagnostics pass, change it
% to 'publication' for the final manuscript results.
%
% Required toolbox: MATLAB Optimization Toolbox (fmincon).
%
% Input (the program searches automatically):
%   trace_driven_inputs_120.csv
%   case_package/outputs/trace_driven_inputs_120.csv
%   case_package/trace_driven_inputs_120.csv
%
% Main outputs (written to trace_case_outputs):
%   Figure_TraceDrivenCase_<mode>.pdf/png
%   Table_TraceCase_StrategyComparison_<mode>.csv
%   Table_TraceCase_Changes_<mode>.csv

clearvars; clc; close all;
rng(20261026,'twister');

runMode = 'publication';  % 'quick' or 'publication'

thisDir = fileparts(mfilename('fullpath'));
inputCandidates = { ...
    fullfile(thisDir,'trace_driven_inputs_120.csv'), ...
    fullfile(thisDir,'case_package','outputs','trace_driven_inputs_120.csv'), ...
    fullfile(thisDir,'case_package','trace_driven_inputs_120.csv')};
inputFile = '';
for iFile = 1:numel(inputCandidates)
    if exist(inputCandidates{iFile},'file') == 2
        inputFile = inputCandidates{iFile};
        break;
    end
end
if isempty(inputFile)
    error(['Trace input file not found. Put trace_driven_inputs_120.csv ' ...
        'in the same folder as Section6_5_TraceDrivenCase.m, or retain ' ...
        'the documented case_package/outputs subfolder.']);
end
outputDir = fullfile(thisDir,'trace_case_outputs');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

p = baseline_parameters(runMode);
p.outputDir = outputDir;
p.inputFile = inputFile;
p = load_trace_workload(p,inputFile);
validate_parameters(p);

if exist('fmincon','file') ~= 2
    error('This program requires fmincon from MATLAB Optimization Toolbox.');
end

tGrid = linspace(0,p.T,p.N+1)';
options = build_solver_options(p);

fprintf('\n===============================================================\n');
fprintf('Trace-driven representative case: %s mode\n',upper(p.runMode));
fprintf('===============================================================\n');
fprintf(['Input=%s\nState nodes=%d, control nodes=%d, starts=%d, ' ...
    'L=0,...,%d\n'],p.inputFile,p.N+1,p.Ncontrol,p.nStarts,p.Lmax);
fprintf('Trace workload: mean=%.8f, min=%.4f, max=%.4f, CV=%.4f\n', ...
    mean(p.yValue),min(p.yValue),max(p.yValue), ...
    std(p.yValue,1)/mean(p.yValue));

%% 1. Constant-workload economic benchmark

pConstant = p;
pConstant.yMode = 'constant';
pConstant.A0 = p.A0Binding;

fprintf('\n--- Case 1: constant workload, soft constraint ---\n');
resConstantSoft = solve_dual_model(pConstant,tGrid,options,false,[]);

%% 2. Trace-driven economic optimum

pTrace = p;
pTrace.yMode = 'trace';
pTrace.A0 = p.A0Binding;

fprintf('\n--- Case 2: trace workload, soft constraint ---\n');
resTraceSoft = solve_dual_model(pTrace,tGrid,options,false,resConstantSoft);

%% 3. Trace-driven hard carbon-account model

fprintf('\n--- Case 3: trace workload, hard constraint ---\n');
if resTraceSoft.minA >= -p.hardAcceptTol
    % Feasible-set nesting: a dynamically safe soft optimum is also optimal
    % for the hard problem, so a redundant hard solve is unnecessary.
    resTraceHard = resTraceSoft;
    hardStatus = 'Nonbinding';
    fprintf(['Trace soft optimum already satisfies A(t)>=0; ' ...
        'the hard solution is identical and compliance cost is zero.\n']);
else
    resTraceHard = solve_dual_model(pTrace,tGrid,options,true,resTraceSoft);
    hardStatus = 'Binding';
    if ~resTraceHard.feasible
        error(['No feasible hard trace solution was obtained. Increase ' ...
            'A0Binding slightly or switch to publication mode.']);
    end
end

% Numerical nesting diagnostic. The soft feasible set contains the hard set.
if resTraceHard.J > resTraceSoft.J+p.valueTol
    warning(['The first trace run returned V^H>V^B. Re-polishing the soft ' ...
        'model once from the hard solution.']);
    resSoftCheck = solve_dual_model(pTrace,tGrid,options,false,resTraceHard);
    if resSoftCheck.J > resTraceSoft.J
        resTraceSoft = resSoftCheck;
    end
end

Ccomp = resTraceSoft.J-resTraceHard.J;
Pcomp = 100*Ccomp/max(abs(resTraceSoft.J),1e-10);
if Ccomp < -p.valueTol
    warning(['Soft-hard value ordering still fails by %.3e. Do not report ' ...
        'quick-mode values; rerun in publication mode.'],-Ccomp);
end
if resTraceHard.minA < -p.hardAcceptTol
    error('Hard trace result violates A(t)>=0 by %.3e.',-resTraceHard.minA);
end

%% 4. Strategy metrics and result tables

mConstant = trapz(tGrid,resConstantSoft.mGrid);
mTraceSoft = trapz(tGrid,resTraceSoft.mGrid);
mTraceHard = trapz(tGrid,resTraceHard.mGrid);

omConstant = sum(resConstantSoft.omega);
omTraceSoft = sum(resTraceSoft.omega);
omTraceHard = sum(resTraceHard.omega);

sigConstant = weighted_upgrade_time( ...
    resConstantSoft.sigma,resConstantSoft.omega);
sigTraceSoft = weighted_upgrade_time( ...
    resTraceSoft.sigma,resTraceSoft.omega);
sigTraceHard = weighted_upgrade_time( ...
    resTraceHard.sigma,resTraceHard.omega);

CaseName = {'Constant';'Trace';'Trace'};
Constraint = {'Soft';'Soft';'Hard'};
RunMode = repmat({p.runMode},3,1);
HardConstraintStatus = {'Not imposed';'Not imposed';hardStatus};
WorkloadMean = [1;mean(p.yValue);mean(p.yValue)];
WorkloadCV = [0;std(p.yValue,1)/mean(p.yValue); ...
                std(p.yValue,1)/mean(p.yValue)];
ObjectiveValue = [resConstantSoft.J;resTraceSoft.J;resTraceHard.J];
MinimumAccount = [resConstantSoft.minA;resTraceSoft.minA;resTraceHard.minA];
TerminalAccount = [resConstantSoft.terminalA; ...
                   resTraceSoft.terminalA;resTraceHard.terminalA];
UpgradeCount = [resConstantSoft.L;resTraceSoft.L;resTraceHard.L];
MaintenanceTotal = [mConstant;mTraceSoft;mTraceHard];
UpgradeIntensityTotal = [omConstant;omTraceSoft;omTraceHard];
WeightedUpgradeTime = [sigConstant;sigTraceSoft;sigTraceHard];
CumulativePrivacyRisk = [resConstantSoft.riskCum; ...
                         resTraceSoft.riskCum;resTraceHard.riskCum];
AverageEfficiency = [resConstantSoft.kappaAvg; ...
                     resTraceSoft.kappaAvg;resTraceHard.kappaAvg];
CumulativeEmissions = [resConstantSoft.carbonCum; ...
                       resTraceSoft.carbonCum;resTraceHard.carbonCum];
ComplianceCost = [NaN;Ccomp;Ccomp];
RelativeComplianceCost = [NaN;Pcomp;Pcomp];

strategyTable = table(CaseName,Constraint,RunMode,HardConstraintStatus, ...
    WorkloadMean,WorkloadCV,ObjectiveValue,MinimumAccount,TerminalAccount, ...
    UpgradeCount,MaintenanceTotal,UpgradeIntensityTotal, ...
    WeightedUpgradeTime,CumulativePrivacyRisk,AverageEfficiency, ...
    CumulativeEmissions,ComplianceCost,RelativeComplianceCost);

comparisonA = build_change_block( ...
    'Trace soft - Constant soft', ...
    resConstantSoft,resTraceSoft,mConstant,mTraceSoft, ...
    omConstant,omTraceSoft,sigConstant,sigTraceSoft);
comparisonB = build_change_block( ...
    'Trace hard - Trace soft', ...
    resTraceSoft,resTraceHard,mTraceSoft,mTraceHard, ...
    omTraceSoft,omTraceHard,sigTraceSoft,sigTraceHard);
changeTable = [comparisonA;comparisonB];

modeTag = lower(p.runMode);
strategyFile = fullfile(outputDir, ...
    sprintf('Table_TraceCase_StrategyComparison_%s.csv',modeTag));
changeFile = fullfile(outputDir, ...
    sprintf('Table_TraceCase_Changes_%s.csv',modeTag));
writetable(strategyTable,strategyFile);
writetable(changeTable,changeFile);

%% 5. Publication-oriented figure

set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',11);

constantColor = [0.45 0.45 0.45];
traceColor    = [0.00 0.32 0.58];
hardColor     = [0.24 0.52 0.30];

% 增大 figure 高度，给顶部 LaTeX 总标题留出空间
fig = figure('Name','Trace-driven representative case', ...
    'Color','w','Position',[60,30,1200,960]);

% 使用 'loose' padding，避免 tiledlayout 总标题被裁切（关键修复）
tl = tiledlayout(fig,2,2,'Padding','loose','TileSpacing','compact');

% (a) Observed workload path.
ax1 = nexttile(tl,[1 2]);
plot(ax1,p.yTime,p.yValue,'-','Color',traceColor,'LineWidth',1.8);
hold(ax1,'on');
yline(ax1,1,':','Color',[0.35 0.35 0.35],'LineWidth',1.2);
hold(ax1,'off');
xlim(ax1,[0 p.T]);
xlabel(ax1,'Time');
ylabel(ax1,'Normalized workload, $y(t)$','Interpreter','latex');
title(ax1,'(a) Trace-derived computing workload');
grid(ax1,'on'); box(ax1,'on'); set(ax1,'LineWidth',1.0);

% (b) Carbon-account paths.
ax2 = nexttile(tl);
plot(ax2,resConstantSoft.pathTime,resConstantSoft.Apath,'--', ...
    'Color',constantColor,'LineWidth',1.6);
hold(ax2,'on');
plot(ax2,resTraceSoft.pathTime,resTraceSoft.Apath,'-', ...
    'Color',traceColor,'LineWidth',1.8);
plot(ax2,resTraceHard.pathTime,resTraceHard.Apath,'-.', ...
    'Color',hardColor,'LineWidth',1.9);
yline(ax2,0,':','Color',[0.25 0.25 0.25],'LineWidth',1.1);
hold(ax2,'off');
xlim(ax2,[0 p.T]);
xlabel(ax2,'Time');
ylabel(ax2,'$A(t)$','Interpreter','latex');
title(ax2,'(b) Dynamic carbon-account paths');
legend(ax2,{'Constant: soft','Trace: soft','Trace: hard'}, ...
    'Location','best','Box','off');
grid(ax2,'on'); box(ax2,'on'); set(ax2,'LineWidth',1.0);

% (c) Trace-governance adjustment after imposing A(t)>=0.
ax3 = nexttile(tl);
yyaxis(ax3,'left');
hmSoft = plot(ax3,tGrid,resTraceSoft.mGrid,'-', ...
    'Color',traceColor,'LineWidth',1.8);
hold(ax3,'on');
hmHard = plot(ax3,tGrid,resTraceHard.mGrid,'-.', ...
    'Color',hardColor,'LineWidth',1.9);
ylabel(ax3,'$m^*(t)$','Interpreter','latex');
ylim(ax3,[0,max(p.mbar,1.05*max([resTraceSoft.mGrid; ...
                                resTraceHard.mGrid]))]);
ax3.YColor = [0 0 0];

yyaxis(ax3,'right');
huSoft = plot_upgrade_stems(ax3,resTraceSoft.sigma, ...
    resTraceSoft.omega,'o',traceColor,false);
huHard = plot_upgrade_stems(ax3,resTraceHard.sigma, ...
    resTraceHard.omega,'s',hardColor,true);
ylabel(ax3,'Upgrade intensity, $\omega_l$','Interpreter','latex');
ylim(ax3,[0 1.05]);
ax3.YColor = [0.30 0.30 0.30];
hold(ax3,'off');
xlim(ax3,[0 p.T]);
xlabel(ax3,'Time');
title(ax3,'(c) Governance response under the hard constraint');
legend(ax3,[hmSoft,hmHard,huSoft,huHard], ...
    {'Maintenance: soft','Maintenance: hard', ...
     'Upgrade: soft','Upgrade: hard'}, ...
    'Location','best','Box','off','FontSize',8);
grid(ax3,'on'); box(ax3,'on'); set(ax3,'LineWidth',1.0);


% 先渲染完成再导出
drawnow;
drawnow;

figurePdf = fullfile(outputDir, ...
    sprintf('Figure_TraceDrivenCase_%s.pdf',modeTag));
figurePng = fullfile(outputDir, ...
    sprintf('Figure_TraceDrivenCase_%s.png',modeTag));

% 导出：显式白底，避免透明区域导致裁剪差异
exportgraphics(fig,figurePdf, ...
    'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,figurePng, ...
    'Resolution',600,'BackgroundColor','white');

% 恢复全局默认字体设置，避免影响后续 figure
set(groot,'defaultAxesFontName','remove');
set(groot,'defaultTextFontName','remove');
set(groot,'defaultAxesFontSize','remove');

%% 6. Console report

fprintf('\n===============================================================\n');
fprintf('Trace-driven comparison summary (%s mode)\n',upper(p.runMode));
fprintf('===============================================================\n');
print_result('Constant / soft',resConstantSoft,mConstant, ...
    omConstant,sigConstant);
print_result('Trace / soft',resTraceSoft,mTraceSoft, ...
    omTraceSoft,sigTraceSoft);
print_result('Trace / hard',resTraceHard,mTraceHard, ...
    omTraceHard,sigTraceHard);
fprintf('Hard-constraint status             = %s\n',hardStatus);
fprintf('Absolute compliance cost Ccomp     = %12.8f\n',Ccomp);
fprintf('Relative compliance cost Pcomp     = %12.6f %%\n',Pcomp);
fprintf('\nFiles written to %s:\n',outputDir);
fprintf('  %s\n',figurePdf);
fprintf('  %s\n',figurePng);
fprintf('  %s\n',strategyFile);
fprintf('  %s\n',changeFile);
if strcmpi(p.runMode,'quick')
    fprintf(['\nQUICK MODE RESULTS ARE DIAGNOSTIC. If feasibility and value ' ...
        'ordering pass, set runMode=''publication'' before reporting.\n']);
end

end


%% ======================================================================
% Parameters and trace input
%% ======================================================================

function p = baseline_parameters(runMode)

p.T  = 10;
p.N  = 120;
p.dt = p.T/p.N;
p.r  = 0.05;

p.runMode = validatestring(runMode,{'quick','publication'});
switch lower(p.runMode)
    case 'quick'
        p.Ncontrol = 10;
        p.nStarts = 1;
        p.maxIterations = 250;
        p.maxFunctionEvaluations = 3e4;
        p.optimalityTolerance = 1e-5;
        p.constraintTolerance = 1e-6;
    case 'publication'
        p.Ncontrol = 16;
        p.nStarts = 3;
        p.maxIterations = 500;
        p.maxFunctionEvaluations = 8e4;
        p.optimalityTolerance = 1e-6;
        p.constraintTolerance = 1e-7;
end
p.stepTolerance = 1e-8;
p.solverDisplay = 'off';

p.yMode = 'constant';
p.yConstant = 1.00;
p.yTime = [];
p.yValue = [];

p.price = 3.00;
p.pE    = 0.60;
p.chi   = 0.35;
p.Q     = 0.50;

p.R0       = 0.40;
p.kappa0   = 0.55;
p.kappaBar = 1.00;
p.A0Safe = 2.20;
p.A0Binding = 0.80;
p.A0 = p.A0Binding;

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

p.Lmax     = 5;
p.sigmaMin = 0.50;
p.sigmaMax = 9.50;
p.minGap   = 0.80;
p.omegaMin = 0.05;

p.bindingTol = 1e-4;
p.hardAcceptTol = 5e-6;
p.valueTol = 1e-6;

end


function p = load_trace_workload(p,inputFile)

if exist(inputFile,'file') ~= 2
    error('Trace input file not found: %s',inputFile);
end

data = readtable(inputFile);
required = {'model_time_start','model_time_end','workload_y'};
if ~all(ismember(required,data.Properties.VariableNames))
    error(['Trace input must contain model_time_start, model_time_end, ' ...
        'and workload_y.']);
end

if height(data) ~= p.N
    error('Trace input must contain exactly %d rows; found %d.', ...
        p.N,height(data));
end

y = data.workload_y(:);
tStart = data.model_time_start(:);
tEnd = data.model_time_end(:);
tMid = 0.5*(tStart+tEnd);

if any(~isfinite(y)) || any(y<=0)
    error('workload_y must contain finite positive values.');
end
if any(~isfinite(tMid)) || any(diff(tMid)<=0)
    error('Trace model times must be finite and strictly increasing.');
end
if abs(mean(y)-1)>1e-6
    error('Trace workload mean must equal one; observed %.10f.',mean(y));
end
if tStart(1)<-1e-9 || tEnd(end)>p.T+1e-9
    error('Trace time interval must lie inside [0,T].');
end

p.yTime = tMid;
p.yValue = y;

end


function options = build_solver_options(p)

options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display',p.solverDisplay, ...
    'MaxIterations',p.maxIterations, ...
    'MaxFunctionEvaluations',p.maxFunctionEvaluations, ...
    'StepTolerance',p.stepTolerance, ...
    'OptimalityTolerance',p.optimalityTolerance, ...
    'ConstraintTolerance',p.constraintTolerance);

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
if ~ismember(lower(p.yMode),{'constant','trace'})
    error('yMode must be constant or trace.');
end
if numel(p.yTime)~=p.N || numel(p.yValue)~=p.N
    error('Trace workload vectors must contain N values.');
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
        [R,K,A] = propagate_state(R,K,A,mMid,t0,dt,p);

        kEff = max(0.5*(K0state+K),1e-6);
        Ravg = 0.5*(R0state+R);
        Kavg = 0.5*(K0state+K);

        yMid = workload_at(tm,p);
        operatingCost = p.pE*yMid/kEff+0.5*p.ay*(yMid/kEff)^2;
        maintenanceCost = p.cm*mMid+0.5*p.am*mMid^2;
        privacyCost = 0.5*p.aR*Ravg^2;
        currentProfit = p.price*yMid-operatingCost-maintenanceCost-privacyCost;

        runningProfit = runningProfit+exp(-p.r*tm)*dt*currentProfit;

        currentCarbon = p.chi*yMid/kEff+p.em*mMid+0.5*p.em2*mMid^2;
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


function [Rnext,Knext,Anext] = propagate_state(R,K,A,m,t0,dt,p)

x = [R;K;A];
k1 = state_rhs(x,m,t0,p);
k2 = state_rhs(x+0.5*dt*k1,m,t0+0.5*dt,p);
k3 = state_rhs(x+0.5*dt*k2,m,t0+0.5*dt,p);
k4 = state_rhs(x+dt*k3,m,t0+dt,p);
xNext = x+dt*(k1+2*k2+2*k3+k4)/6;

Rnext = max(xNext(1),0);
Knext = min(max(xNext(2),1e-6),p.kappaBar);
Anext = xNext(3);

end


function dx = state_rhs(x,m,t,p)

R = max(x(1),0);
K = min(max(x(2),1e-6),p.kappaBar);
y = workload_at(t,p);

dR = p.beta*y-(p.deltaR+p.etaR*m)*R;
dK = p.etaK*m*(p.kappaBar-K)-p.deltaK*K;
dA = p.Q-p.chi*y/K-p.em*m-0.5*p.em2*m^2;

dx = [dR;dK;dA];

end


function y = workload_at(t,p)

switch lower(p.yMode)
    case 'constant'
        y = p.yConstant;
    case 'trace'
        tq = min(max(t,p.yTime(1)),p.yTime(end));
        y = interp1(p.yTime,p.yValue,tq,'linear');
    otherwise
        error('Unknown workload mode: %s',p.yMode);
end

end


%% ======================================================================
% Reporting helpers
%% ======================================================================

function changeTable = build_change_block(label,res0,res1, ...
    m0,m1,omega0,omega1,sigma0,sigma1)

Comparison = repmat({label},8,1);
Indicator = {'Objective value';'Minimum account';'Maintenance'; ...
    'Upgrade intensity';'Weighted upgrade time';'Privacy risk'; ...
    'Average efficiency';'Cumulative emissions'};
BaseValue = [res0.J;res0.minA;m0;omega0;sigma0;res0.riskCum; ...
             res0.kappaAvg;res0.carbonCum];
ComparisonValue = [res1.J;res1.minA;m1;omega1;sigma1;res1.riskCum; ...
                   res1.kappaAvg;res1.carbonCum];
AbsoluteChange = ComparisonValue-BaseValue;
RelativeChangePercent = nan(size(BaseValue));
for i = 1:numel(BaseValue)
    RelativeChangePercent(i) = relative_change( ...
        ComparisonValue(i),BaseValue(i));
end

changeTable = table(Comparison,Indicator,BaseValue,ComparisonValue, ...
    AbsoluteChange,RelativeChangePercent);

end


function h = plot_upgrade_stems(ax,sigma,omega,marker,color,filled)

if isempty(sigma)
    h = plot(ax,nan,nan,marker,'Color',color,'LineStyle','none');
    return;
end

h = stem(ax,sigma,omega,marker,'Color',color,'LineWidth',1.1, ...
    'MarkerSize',5);
if filled
    h.MarkerFaceColor = color;
else
    h.MarkerFaceColor = 'w';
end

end


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
fprintf('%-18s cumulative emissions=%10.6f\n', ...
    name,res.carbonCum);

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