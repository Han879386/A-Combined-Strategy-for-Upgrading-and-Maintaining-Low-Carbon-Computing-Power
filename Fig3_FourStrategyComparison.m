function Figure3_FourStrategyComparison
%% ========================================================================
% Endogenous maintenance, upgrade timing, intensity, and upgrade count
%
% Active upgrades satisfy omega in [omegaMin,1]. The no-upgrade decision is
% represented by L=0. The active upgrade count is selected from
% L=0,1,...,Lfeas, where Lfeas is implied by the finite upgrade window and
% the minimum spacing condition.
%
% Standalone Figure 3: performance comparison among four strategies
%
% A fixed setup cost K0 is charged whenever an active major upgrade occurs.
% This cost affects whether and when to upgrade, but it does not enter the
% interior first-order condition for upgrade intensity.
%% ========================================================================

clear; clc; close all;

p = baseline_parameters();
tGrid = linspace(0,p.T,p.N+1)';

if exist('fmincon','file') ~= 2
    error('Optimization Toolbox is required because this code uses fmincon.');
end

options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','off', ...
    'MaxIterations',1800, ...
    'MaxFunctionEvaluations',1e6, ...
    'StepTolerance',1e-9, ...
    'OptimalityTolerance',1e-7, ...
    'ConstraintTolerance',1e-8);

%% Natural count bound implied by timing restrictions
Lfeas = floor((p.sigmaMax-p.sigmaMin+1e-12)/p.minGap)+1;
if isempty(p.countSearchCap)
    Lsearch = Lfeas;
else
    Lsearch = min(p.countSearchCap,Lfeas);
end

fprintf('\nNatural feasible upgrade count: %d\n',Lfeas);
fprintf('Numerically searched counts: 0,...,%d\n',Lsearch);
fprintf('Minimum effective upgrade intensity: %.4f\n',p.omegaMin);
fprintf('Fixed setup cost per active upgrade K0: %.4f\n',p.K0);

if Lsearch < Lfeas
    warning('The count search is truncated. Set p.countSearchCap=[] for the final full search.');
end

%% Solve four strategies
res_none = simulate_policy(zeros(p.N+1,1),[],[],p,tGrid);
res_monly = solve_fixed_count('maintenance',0,p,tGrid,options);
res_uonly = solve_endogenous_count('upgrade',res_none,Lsearch,p,tGrid,options);
res_dual = solve_endogenous_count('dual',res_monly,Lsearch,p,tGrid,options);

%% Figures
plot_strategy_comparison(res_none,res_monly,res_uonly,res_dual);

%% Results
fprintf('\n================ Results ================\n');
print_result('No governance',res_none);
print_result('Maintenance only',res_monly);
print_result('Upgrade only',res_uonly);
print_result('Dual policy',res_dual);

fprintf('\nSelected dual-policy count L* = %d\n',res_dual.Lstar);
fprintf('Upgrade times: '); fprintf('%.4f ',res_dual.sigma); fprintf('\n');
fprintf('Upgrade intensities: '); fprintf('%.4f ',res_dual.omega); fprintf('\n');
fprintf('Minimum carbon-account balance = %.6f\n',min(res_dual.Apath));

if res_dual.Lstar == Lsearch
    warning(['The selected count reaches the largest count searched. ' ...
        'Search the full natural range and test minGap and omegaMin.']);
end

if min(res_dual.Apath) < -1e-6
    warning(['The economic optimum produces a negative carbon account at some dates. ' ...
        'This is permitted by the current soft-penalty model and should not be ' ...
        'removed by ad hoc parameter tuning.']);
end

end

%% ========================================================================
% Parameters
%% ========================================================================
function p = baseline_parameters()

p.T = 10;
p.N = 50;                 % follows the current manuscript table
p.r = 0.05;

p.y = 1.00;
p.price = 3.00;
p.pE = 0.60;
p.chi = 0.35;
p.Q = 0.50;
p.deltaA = 0.02;

p.R0 = 0.40;
p.kappa0 = 0.55;
p.A0 = 2.20;
p.kappaBar = 1.00;

p.beta = 0.25;
p.deltaR = 0.05;
p.etaR = 0.45;

p.etaK = 0.35;
p.deltaK = 0.12;

p.cm = 0.08;
p.am = 0.45;
p.em = 0.03;
p.em2 = 0.10;

p.aR = 1.20;
p.ay = 0.45;

p.K0 = 0.05;              % fixed setup cost for each active major upgrade
p.K1 = 0.25;
p.K2 = 0.70;
p.xi1 = 0.14;
p.xi2 = 0.35;
p.gammaR = 0.65;
p.gammaK = 0.65;          % follows the current manuscript table

p.vA = 0.50;
p.phiA = 8.00;
p.vK = 0.80;
p.phiR = 1.50;
p.mbar = 1.00;

p.sigmaMin = 0.50;
p.sigmaMax = 9.50;
p.minGap = 0.80;
p.omegaMin = 0.03;

% [] searches the complete count range implied by timing restrictions.
% For preliminary debugging only, this may be set to a smaller integer.
p.countSearchCap = [];
p.numStarts = 4;

% Counts whose objective is within this tolerance of the maximum are
% reported as a near-optimal count set. This is useful when J_L is flat.
p.nearOptimalRelTol = 1.0e-3;   % 0.10% of |J_max|
p.nearOptimalAbsTol = 1.0e-3;   % absolute lower bound

end

%% ========================================================================
% Endogenous count selection
%% ========================================================================
function bestRes = solve_endogenous_count(strategy,baseRes,Lsearch,p,tGrid,options)

profileL = (0:Lsearch)';
profileJ = nan(Lsearch+1,1);
profileExit = nan(Lsearch+1,1);

profileJ(1) = baseRes.J;
profileExit(1) = 1;
bestRes = baseRes;
bestJ = baseRes.J;

fprintf('\nSearching active upgrade count for %s\n',strategy);

for L = 1:Lsearch
    fprintf('  L=%d ... ',L);
    resL = solve_fixed_count(strategy,L,p,tGrid,options);
    profileJ(L+1) = resL.J;
    profileExit(L+1) = resL.exitflag;
    fprintf('J=%.8f, exitflag=%d\n',resL.J,resL.exitflag);

    if resL.J > bestJ
        bestJ = resL.J;
        bestRes = resL;
    end
end

bestRes.Lstar = numel(bestRes.omega);
bestRes.countProfile.L = profileL;
bestRes.countProfile.J = profileJ;
bestRes.countProfile.exitflag = profileExit;

end

%% ========================================================================
% Fixed active-upgrade count problem
%% ========================================================================
function res = solve_fixed_count(strategy,L,p,tGrid,options)

N = p.N;

switch lower(strategy)
    case 'maintenance'
        nvar = N+1;
        lb = zeros(nvar,1);
        ub = p.mbar*ones(nvar,1);
        Aineq = [];
        bineq = [];
        levels = [0.10,0.35,0.65,0.90];
        z0_set = zeros(nvar,numel(levels));
        for k=1:numel(levels)
            z0_set(:,k) = levels(k)*ones(nvar,1);
        end

    case 'upgrade'
        nvar = 2*L;
        lb = [p.sigmaMin*ones(L,1);p.omegaMin*ones(L,1)];
        ub = [p.sigmaMax*ones(L,1);ones(L,1)];
        [Aineq,bineq] = spacing_constraints(1,L,nvar,p.minGap);
        z0_set = make_initial_points('upgrade',L,p);

    case 'dual'
        nvar = (N+1)+2*L;
        lb = [zeros(N+1,1);p.sigmaMin*ones(L,1);p.omegaMin*ones(L,1)];
        ub = [p.mbar*ones(N+1,1);p.sigmaMax*ones(L,1);ones(L,1)];
        sigmaStart = N+2;
        [Aineq,bineq] = spacing_constraints(sigmaStart,L,nvar,p.minGap);
        z0_set = make_initial_points('dual',L,p);

    otherwise
        error('Unknown strategy.');
end

bestObj = Inf;
bestZ = [];
bestExitflag = -Inf;
bestOutput = struct();

nStarts = min(size(z0_set,2),p.numStarts);
for k=1:nStarts
    z0 = z0_set(:,k);
    objFun = @(z)objective_wrapper(z,strategy,L,p,tGrid);
    [zSol,objVal,exitflag,output] = fmincon( ...
        objFun,z0,Aineq,bineq,[],[],lb,ub,[],options);

    if objVal < bestObj
        bestObj = objVal;
        bestZ = zSol;
        bestExitflag = exitflag;
        bestOutput = output;
    end
end

if isempty(bestZ)
    error('No solution returned for %s with L=%d.',strategy,L);
end

[mGrid,sigma,omega] = unpack_controls(bestZ,strategy,L,p);
res = simulate_policy(mGrid,sigma,omega,p,tGrid);
res.exitflag = bestExitflag;
res.output = bestOutput;
res.Lstar = L;
res.localDiagnostics = local_event_checks(res,p,tGrid);

end

function [Aineq,bineq] = spacing_constraints(sigmaStart,L,nvar,minGap)
if L<=1
    Aineq=[]; bineq=[]; return;
end
Aineq=zeros(L-1,nvar);
bineq=-minGap*ones(L-1,1);
for i=1:L-1
    Aineq(i,sigmaStart+i-1)=1;
    Aineq(i,sigmaStart+i)=-1;
end
end

function z0_set = make_initial_points(strategy,L,p)
S = feasible_sigma_starts(L,p);
omegaLevels=[0.25,0.45,0.65,0.85];
mLevels=[0.20,0.40,0.65,0.85];
nStarts=min([size(S,2),numel(omegaLevels),p.numStarts]);

switch lower(strategy)
    case 'upgrade'
        z0_set=zeros(2*L,nStarts);
        for k=1:nStarts
            z0_set(:,k)=[S(:,k);max(p.omegaMin,omegaLevels(k))*ones(L,1)];
        end
    case 'dual'
        z0_set=zeros((p.N+1)+2*L,nStarts);
        for k=1:nStarts
            z0_set(:,k)=[mLevels(k)*ones(p.N+1,1);S(:,k); ...
                max(p.omegaMin,omegaLevels(k))*ones(L,1)];
        end
end
end

function S = feasible_sigma_starts(L,p)
if L==1
    S=[p.sigmaMin,0.5*(p.sigmaMin+p.sigmaMax),p.sigmaMax, ...
        p.sigmaMin+0.35*(p.sigmaMax-p.sigmaMin)];
    return;
end
early=p.sigmaMin+(0:L-1)'*p.minGap;
late=p.sigmaMax-(L-1:-1:0)'*p.minGap;
even=linspace(p.sigmaMin,p.sigmaMax,L)';
slack=p.sigmaMax-p.sigmaMin-(L-1)*p.minGap;
centerStart=p.sigmaMin+0.5*max(slack,0);
centered=centerStart+(0:L-1)'*p.minGap;
S=[early,even,late,centered];
end

%% ========================================================================
% Objective and unpacking
%% ========================================================================
function obj = objective_wrapper(z,strategy,L,p,tGrid)
[mGrid,sigma,omega]=unpack_controls(z,strategy,L,p);
res=simulate_policy(mGrid,sigma,omega,p,tGrid);
obj=-res.J;
end

function [mGrid,sigma,omega]=unpack_controls(z,strategy,L,p)
N=p.N;
switch lower(strategy)
    case 'maintenance'
        mGrid=z(1:N+1); sigma=[]; omega=[];
    case 'upgrade'
        mGrid=zeros(N+1,1);
        sigma=z(1:L); omega=z(L+1:2*L);
    case 'dual'
        mGrid=z(1:N+1);
        sigma=z(N+2:N+1+L);
        omega=z(N+2+L:N+1+2*L);
end
end

%% ========================================================================
% Simulation
%% ========================================================================
function res=simulate_policy(mGrid,sigma,omega,p,tGrid)

sigma=sigma(:); omega=omega(:);
if numel(sigma)~=numel(omega)
    error('sigma and omega must have the same length.');
end
if ~isempty(omega) && any(omega<p.omegaMin-1e-9)
    error('Every active upgrade must satisfy omega>=omegaMin.');
end

allTimes=sort([tGrid(:);sigma(:)]);
timeList=allTimes([true;diff(allTimes)>1e-10]);

R=p.R0; K=p.kappa0; A=p.A0;
pathTime=0; Rpath=R; Kpath=K; Apath=A;
runningProfit=0; runningCarbon=0; riskCum=0; kappaCum=0;
upgradeCostPV=0; upgradeFixedCostPV=0; upgradeVariableCostPV=0;
upgradeCarbon=0;
eventApplied=false(length(sigma),1);

for q=1:length(timeList)-1
    t0=timeList(q); t1=timeList(q+1); dt=t1-t0;

    if dt>0
        tm=0.5*(t0+t1);
        mMid=interp1(tGrid,mGrid,tm,'linear','extrap');
        mMid=min(max(mMid,0),p.mbar);

        Rprev=R; Kprev=K; Aprev=A;
        [R,K,A]=propagate_state(R,K,A,mMid,dt,p);

        kEff=max(0.5*(Kprev+K),1e-8);
        Ravg=0.5*(Rprev+R);
        Kavg=0.5*(Kprev+K);

        operatingCost=p.pE*p.y/kEff+0.5*p.ay*(p.y/kEff)^2;
        maintenanceCost=p.cm*mMid+0.5*p.am*mMid^2;
        privacyCost=0.5*p.aR*Ravg^2;
        currentProfit=p.price*p.y-operatingCost-maintenanceCost-privacyCost;
        runningProfit=runningProfit+exp(-p.r*tm)*dt*currentProfit;

        currentCarbon=p.chi*p.y/kEff+p.em*mMid+0.5*p.em2*mMid^2;
        runningCarbon=runningCarbon+dt*currentCarbon;
        riskCum=riskCum+dt*Ravg;
        kappaCum=kappaCum+dt*Kavg;

        pathTime(end+1,1)=t1;
        Rpath(end+1,1)=R;
        Kpath(end+1,1)=K;
        Apath(end+1,1)=A;
    end

    idxEvents=find(abs(sigma-t1)<1e-8 & ~eventApplied);
    for jj=1:length(idxEvents)
        l=idxEvents(jj); w=omega(l);
        Kminus=K; Rminus=R; Aminus=A;
        K=Kminus+p.gammaK*w*(p.kappaBar-Kminus);
        R=(1-p.gammaR*w)*Rminus;
        A=Aminus-p.xi1*w-0.5*p.xi2*w^2;

        discountFactor=exp(-p.r*sigma(l));
        fixedCostPV=discountFactor*p.K0;
        variableCostPV=discountFactor*(p.K1*w+0.5*p.K2*w^2);

        upgradeFixedCostPV=upgradeFixedCostPV+fixedCostPV;
        upgradeVariableCostPV=upgradeVariableCostPV+variableCostPV;
        upgradeCostPV=upgradeCostPV+fixedCostPV+variableCostPV;
        upgradeCarbon=upgradeCarbon+p.xi1*w+0.5*p.xi2*w^2;

        pathTime(end+1,1)=t1;
        Rpath(end+1,1)=R;
        Kpath(end+1,1)=K;
        Apath(end+1,1)=A;
        eventApplied(l)=true;
    end
end

Aplus=max(A,0); Aminus=max(-A,0);
terminalValue=p.vA*Aplus-0.5*p.phiA*Aminus^2+ ...
    p.vK*K-0.5*p.phiR*R^2;
J=runningProfit-upgradeCostPV+exp(-p.r*p.T)*terminalValue;

res.mGrid=mGrid;
res.sigma=sigma;
res.omega=omega;
res.pathTime=pathTime;
res.Rpath=Rpath;
res.Kpath=Kpath;
res.Apath=Apath;
res.J=J;
res.riskCum=riskCum;
res.carbonCum=runningCarbon+upgradeCarbon;
res.kappaAvg=kappaCum/p.T;
res.upgradeCostPV=upgradeCostPV;
res.upgradeFixedCostPV=upgradeFixedCostPV;
res.upgradeVariableCostPV=upgradeVariableCostPV;
res.exitflag=1;
res.Lstar=numel(omega);
end

function [Rnext,Knext,Anext]=propagate_state(R,K,A,m,dt,p)
x=[R;K;A]; f=@(xv)state_rhs(xv,m,p);
k1=f(x); k2=f(x+0.5*dt*k1); k3=f(x+0.5*dt*k2); k4=f(x+dt*k3);
xNext=x+dt*(k1+2*k2+2*k3+k4)/6;
Rnext=max(xNext(1),0);
Knext=min(max(xNext(2),1e-8),p.kappaBar);
Anext=xNext(3);
end

function dx=state_rhs(x,m,p)
R=max(x(1),0); K=min(max(x(2),1e-8),p.kappaBar); A=x(3);
dR=p.beta*p.y-(p.deltaR+p.etaR*m)*R;
dK=p.etaK*m*(p.kappaBar-K)-p.deltaK*K;
dA=p.Q-p.chi*p.y/K-p.em*m-0.5*p.em2*m^2-p.deltaA*A;
dx=[dR;dK;dA];
end

%% ========================================================================
% Local perturbation diagnostics
%% ========================================================================
function diag=local_event_checks(res,p,tGrid)
L=numel(res.omega);
diag.timeMinus=nan(L,1); diag.timePlus=nan(L,1);
diag.intMinus=nan(L,1); diag.intPlus=nan(L,1);
if L==0, return; end

epsT=min(0.02,0.1*p.minGap); epsW=0.01;
for l=1:L
    sigma0=res.sigma; omega0=res.omega;
    lowerT=p.sigmaMin; upperT=p.sigmaMax;
    if l>1, lowerT=max(lowerT,sigma0(l-1)+p.minGap); end
    if l<L, upperT=min(upperT,sigma0(l+1)-p.minGap); end

    if sigma0(l)-epsT>=lowerT
        s=sigma0; s(l)=sigma0(l)-epsT;
        r=simulate_policy(res.mGrid,s,omega0,p,tGrid);
        diag.timeMinus(l)=r.J-res.J;
    end
    if sigma0(l)+epsT<=upperT
        s=sigma0; s(l)=sigma0(l)+epsT;
        r=simulate_policy(res.mGrid,s,omega0,p,tGrid);
        diag.timePlus(l)=r.J-res.J;
    end
    if omega0(l)-epsW>=p.omegaMin
        w=omega0; w(l)=omega0(l)-epsW;
        r=simulate_policy(res.mGrid,sigma0,w,p,tGrid);
        diag.intMinus(l)=r.J-res.J;
    end
    if omega0(l)+epsW<=1
        w=omega0; w(l)=omega0(l)+epsW;
        r=simulate_policy(res.mGrid,sigma0,w,p,tGrid);
        diag.intPlus(l)=r.J-res.J;
    end
end
end

%% ========================================================================
% Figures
%% ========================================================================
function plot_strategy_comparison(res_none,res_monly,res_uonly,res_dual)
strategyNames={'No governance','Maintenance','Upgrade','Dual'};
J=[res_none.J,res_monly.J,res_uonly.J,res_dual.J];
Kappa=[res_none.kappaAvg,res_monly.kappaAvg,res_uonly.kappaAvg,res_dual.kappaAvg];
Carbon=[res_none.carbonCum,res_monly.carbonCum,res_uonly.carbonCum,res_dual.carbonCum];
Risk=[res_none.riskCum,res_monly.riskCum,res_uonly.riskCum,res_dual.riskCum];

Y={J,Kappa,Carbon,Risk};
Titles={'Discounted profit','Average efficiency', ...
    'Carbon emissions','Privacy risk'};
Labels={'$J$','$\bar{\kappa}$','$C_T$', ...
    '$\int_0^T R(t)\,dt$'};

figure('Name','Figure 3: Strategy comparison','Color','w', ...
    'Position',[80,180,1450,390],'ToolBar','none');
tiledlayout(1,4,'Padding','compact','TileSpacing','compact');

for ii=1:4
    ax=nexttile;
    vals=Y{ii};

    bar(vals,'FaceColor','w','EdgeColor','k','LineWidth',1.2);
    set(ax,'XTick',1:4,'XTickLabel',strategyNames, ...
        'XTickLabelRotation',28);
    ylabel(Labels{ii},'Interpreter','latex');
    title(Titles{ii},'FontWeight','normal');
    grid on; box on;
    set(ax,'LineWidth',1.0,'FontSize',11);

    % Add sufficient headroom so numerical labels do not touch the titles.
    if ii==1
        dataRange=max(vals)-min(vals);
        if dataRange<=0, dataRange=1; end
        lower=min(0,min(vals)-0.08*dataRange);
        upper=max(0,max(vals)+0.16*dataRange);
    else
        upper=max(vals)*1.18;
        lower=0;
    end
    ylim(ax,[lower upper]);

    yl=ylim(ax);
    span=max(yl(2)-yl(1),1);
    offset=0.025*span;

    for j=1:length(vals)
        if vals(j)>=0
            text(ax,j,vals(j)+offset,sprintf('%.3f',vals(j)), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom','FontSize',9);
        else
            text(ax,j,vals(j)-offset,sprintf('%.3f',vals(j)), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','top','FontSize',9);
        end
    end

    clean_axes(ax);
end
end

function clean_axes(ax)
% Remove interactive axes toolbar from screenshots and exported figures.
try
    axtoolbar(ax,{});
catch
end
try
    disableDefaultInteractivity(ax);
catch
end
end

function txt=format_integer_set(values)
values=values(:)';
if isempty(values)
    txt='none';
elseif numel(values)==1
    txt=sprintf('%d',values);
else
    txt=sprintf('%d,',values);
    txt=txt(1:end-1);
end
end

%% ========================================================================
% Reporting
%% ========================================================================
function print_result(name,res)
fprintf('%-20s J=%10.5f Risk=%10.5f Carbon=%10.5f AvgKappa=%10.5f L=%d\n', ...
    name,res.J,res.riskCum,res.carbonCum,res.kappaAvg,numel(res.omega));
fprintf('%-20s Upgrade cost PV: total=%9.5f fixed=%9.5f variable=%9.5f\n', ...
    name,res.upgradeCostPV,res.upgradeFixedCostPV,res.upgradeVariableCostPV);
if isempty(res.sigma)
    fprintf('%-20s sigma=(none)\n',name);
    fprintf('%-20s omega=(none)\n\n',name);
else
    fprintf('%-20s sigma=',name); fprintf('%7.4f ',res.sigma); fprintf('\n');
    fprintf('%-20s omega=',name); fprintf('%7.4f ',res.omega); fprintf('\n\n');
end
end

