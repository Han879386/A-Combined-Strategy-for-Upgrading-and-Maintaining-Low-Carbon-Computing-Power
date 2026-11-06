%% Figure 7: Governance and carbon-account risk under stochastic carbon intensity
% TRUE MATLAB SCRIPT: no user-defined or local functions; no fmincon.
%
% Save this file as Figure7_StochasticCarbonIntensity_Final.m and click Run.
%
% Scientific purpose
% ------------------
% This script implements the stochastic-carbon-intensity extension using a
% scenario-average, OPEN-LOOP policy.  The same maintenance path and major-
% upgrade decisions are used in every scenario at a given volatility level.
% Therefore, changes in maintenance should be interpreted as ex-ante changes
% in the time-divisible governance plan, not as real-time adaptive control.
%
% Carbon-intensity process
%   d chi = muChi*(chiBar-chi)*dt + sigmaChi*dW,
% followed by symmetric winsorization around chiBar.  The symmetric lower
% and upper bounds, together with antithetic shocks, preserve the scenario
% mean at chiBar and isolate volatility from a changing mean intensity.
%
% Final baseline carbon-account model
%   dA/dt = Q - chi(t)*y/kappa(t) - em*m(t) - 0.5*em2*m(t)^2,
%   A(sigma_l+) = A(sigma_l-) - xi1*omega_l - 0.5*xi2*omega_l^2.
%
% Important: the final baseline contains NO natural carbon-account decay
% term (-deltaA*A). Temporary deficits are permitted and penalized only
% through the terminal carbon-account value.
%
% Numerical method
% ----------------
% To keep the extension reproducible and fast, the code searches a common
% policy library rather than repeatedly calling fmincon. Maintenance is a
% piecewise-linear path with four nodes. Each candidate contains 0--4 major
% upgrades whose grid-aligned times satisfy the minimum-spacing condition.
% Risk and efficiency dynamics are deterministic conditional on a policy;
% only the carbon account changes across carbon-intensity scenarios.
%
% Main figure
%   (a) total major-upgrade intensity, Omega*;
%   (b) cumulative maintenance effort;
%   (c) out-of-sample expected objective and 95% confidence interval;
%   (d) dynamic, terminal, and temporary carbon-deficit probabilities.
%
% Additional saved diagnostics
%   effective upgrade count, intensity-weighted upgrade time,
%   expected maximum deficit severity, and portfolio carbon-recovery
%   probability relative to a no-upgrade branch with the SAME maintenance.

clearvars;
clc;
close all;

tic;
rng(20260911,'twister');

%% 1. Structural parameters: consistent with the final benchmark

T = 10.0;
N = 40;                         % fast extension grid
dt = T/N;
tGrid = linspace(0,T,N+1);
tMid = 0.5*(tGrid(1:end-1)+tGrid(2:end));

r = 0.05;
y = 1.00;
price = 3.00;
pE = 0.60;

chiBar = 0.35;
Q = 0.50;

R0 = 0.40;
kappa0 = 0.55;
A0 = 2.20;
kappaBar = 1.00;

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

K0 = 0.08;
K1 = 0.25;
K2 = 0.70;

xi1 = 0.14;
xi2 = 0.35;

gammaR = 0.65;
gammaK = 0.70;

vA = 0.50;
phiA = 8.00;
vK = 0.80;
phiR = 1.50;

mbar = 1.00;

Lmax = 4;
sigmaMin = 0.50;
sigmaMax = 9.50;
minGap = 0.80;
omegaMin = 0.05;

%% 2. Stochastic and numerical settings

muChi = 0.60;
chiLower = 0.01;
chiUpper = 2*chiBar-chiLower;    % 0.69: symmetric around chiBar=0.35
sigmaChiGrid = [0.00 0.05 0.10 0.15 0.20 0.30 0.40 0.50]';
nSigma = numel(sigmaChiGrid);

% Mopt follows the original stochastic extension. Meval is increased to
% provide more stable out-of-sample probabilities and confidence intervals.
Mopt = 35;
Meval = 1000;

% Policy-library size. Raise nRandomPolicy for a publication replication
% check; the current value is designed to finish quickly on ordinary PCs.
nRandomPolicy = 2200;
nMaintenanceNodes = 4;
maintenanceNodeTime = linspace(0,T,nMaintenanceNodes);

% Output settings.
exportFigures = true;
figurePdf = 'Figure7_StochasticCarbonIntensity_Final.pdf';
figurePng = 'Figure7_StochasticCarbonIntensity_Final.png';
resultCsv = 'Figure7_StochasticCarbonIntensity_Results.csv';
resultMat = 'Figure7_StochasticCarbonIntensity_Results.mat';

%% 3. Construct a common open-loop policy library

stepMin = max(1,ceil(sigmaMin/dt));
stepMax = min(N,floor(sigmaMax/dt));
gapSteps = ceil(minGap/dt);

% Random smooth maintenance profiles.
mNodeLibrary = mbar*rand(nRandomPolicy,nMaintenanceNodes);
mNodeLibrary(:,2) = 0.25*mNodeLibrary(:,1) ...
                  + 0.50*mNodeLibrary(:,2) ...
                  + 0.25*mNodeLibrary(:,3);
mNodeLibrary(:,3) = 0.25*mNodeLibrary(:,2) ...
                  + 0.50*mNodeLibrary(:,3) ...
                  + 0.25*mNodeLibrary(:,4);

eventStepLibrary = zeros(nRandomPolicy,Lmax);
omegaLibrary = zeros(nRandomPolicy,Lmax);

for ic = 1:nRandomPolicy
    L = randi([0 Lmax]);

    if L > 0
        slackSteps = stepMax-stepMin-(L-1)*gapSteps;
        randomSlack = sort(floor((slackSteps+1)*rand(1,L)));
        eventStepLibrary(ic,1:L) = stepMin ...
            + randomSlack + (0:L-1)*gapSteps;
        omegaLibrary(ic,1:L) = omegaMin ...
            +(1-omegaMin)*rand(1,L);
    end
end

% Add structured policies so that boundary and constant-control cases are
% always represented rather than left to random sampling.
mLevelSet = 0:0.20:mbar;
omegaLevelSet = [0.20 0.40 0.60 0.80 1.00];

for im = 1:numel(mLevelSet)
    mLevel = mLevelSet(im);

    % No-upgrade candidate.
    mNodeLibrary(end+1,:) = mLevel*ones(1,nMaintenanceNodes);
    eventStepLibrary(end+1,:) = zeros(1,Lmax);
    omegaLibrary(end+1,:) = zeros(1,Lmax);

    for L = 1:Lmax
        earlySteps = stepMin+(0:L-1)*gapSteps;
        spreadSteps = round(linspace(stepMin,stepMax,L));
        lateSteps = stepMax-(L-1:-1:0)*gapSteps;
        timingPatterns = [earlySteps;spreadSteps;lateSteps];

        for ip = 1:size(timingPatterns,1)
            for iw = 1:numel(omegaLevelSet)
                newSteps = zeros(1,Lmax);
                newOmega = zeros(1,Lmax);
                newSteps(1:L) = timingPatterns(ip,:);
                newOmega(1:L) = omegaLevelSet(iw);

                mNodeLibrary(end+1,:) = ...
                    mLevel*ones(1,nMaintenanceNodes);
                eventStepLibrary(end+1,:) = newSteps;
                omegaLibrary(end+1,:) = newOmega;
            end
        end
    end
end

% For each upgrade portfolio, create lower/base/higher maintenance variants.
% This prevents the random library from mechanically tying one upgrade plan
% to only one maintenance path and allows the governance mix to adjust along
% the two margins separately when carbon-intensity volatility changes.
baseMNodeLibrary = mNodeLibrary;
baseEventStepLibrary = eventStepLibrary;
baseOmegaLibrary = omegaLibrary;
nBasePolicy = size(baseMNodeLibrary,1);
maintenanceShiftSet = [-0.15 0.00 0.15];

mNodeLibrary = zeros(nBasePolicy*numel(maintenanceShiftSet), ...
    nMaintenanceNodes);
eventStepLibrary = zeros(nBasePolicy*numel(maintenanceShiftSet),Lmax);
omegaLibrary = zeros(nBasePolicy*numel(maintenanceShiftSet),Lmax);

for ishift = 1:numel(maintenanceShiftSet)
    rowIndex = (ishift-1)*nBasePolicy+(1:nBasePolicy);
    mNodeLibrary(rowIndex,:) = min(max(baseMNodeLibrary ...
        +maintenanceShiftSet(ishift),0),mbar);
    eventStepLibrary(rowIndex,:) = baseEventStepLibrary;
    omegaLibrary(rowIndex,:) = baseOmegaLibrary;
end

nPolicy = size(mNodeLibrary,1);

% Linear interpolation basis for the four maintenance nodes.
maintenanceBasis = zeros(nMaintenanceNodes,N+1);
for j = 1:nMaintenanceNodes
    unitNode = zeros(1,nMaintenanceNodes);
    unitNode(j) = 1;
    maintenanceBasis(j,:) = interp1(maintenanceNodeTime,unitNode, ...
        tGrid,'linear');
end

mLibrary = mNodeLibrary*maintenanceBasis;
mLibrary = min(max(mLibrary,0),mbar);

%% 4. Precompute deterministic state paths and policy values

R = R0*ones(nPolicy,1);
kappa = kappa0*ones(nPolicy,1);

energyCoefficient = zeros(nPolicy,N);       % y/kappa over each interval
maintenanceEmissionRate = zeros(nPolicy,N);
upgradeCarbonStep = zeros(nPolicy,N);

runningValue = zeros(nPolicy,1);
upgradeCostPV = zeros(nPolicy,1);

for n = 1:N
    mMid = 0.5*(mLibrary(:,n)+mLibrary(:,n+1));

    Rold = R;
    kappaOld = kappa;

    % Continuous dynamics over [t_n,t_{n+1}); event is applied at t_{n+1}.
    Rcontinuous = Rold+dt*(beta*y-(deltaR+etaR*mMid).*Rold);
    kappaContinuous = kappaOld+dt*( ...
        etaK*mMid.*(kappaBar-kappaOld)-deltaK*kappaOld);

    Rcontinuous = max(Rcontinuous,0);
    kappaContinuous = min(max(kappaContinuous,1e-8),kappaBar);

    Raverage = 0.5*(Rold+Rcontinuous);
    kappaAverage = max(0.5*(kappaOld+kappaContinuous),1e-8);

    energyCoefficient(:,n) = y./kappaAverage;
    maintenanceEmissionRate(:,n) = ...
        em*mMid+0.5*em2*mMid.^2;

    operatingCost = pE*y./kappaAverage ...
        +0.5*ay*(y./kappaAverage).^2;
    maintenanceCost = cm*mMid+0.5*am*mMid.^2;
    privacyCost = 0.5*aR*Raverage.^2;

    currentProfit = price*y-operatingCost ...
        -maintenanceCost-privacyCost;
    runningValue = runningValue ...
        +exp(-r*tMid(n))*dt*currentProfit;

    eventMask = eventStepLibrary==n;
    omegaNow = sum(omegaLibrary.*eventMask,2);
    activeNow = omegaNow>0;

    upgradeCarbonStep(:,n) = ...
        xi1*omegaNow+0.5*xi2*omegaNow.^2;

    upgradeCostPV = upgradeCostPV ...
        +exp(-r*tGrid(n+1)).*( ...
        K0*double(activeNow)+K1*omegaNow+0.5*K2*omegaNow.^2);

    % State jumps from the major upgrade.
    R = (1-gammaR*omegaNow).*Rcontinuous;
    kappa = kappaContinuous ...
        +gammaK*omegaNow.*(kappaBar-kappaContinuous);

    R = max(R,0);
    kappa = min(max(kappa,1e-8),kappaBar);
end

discountTerminal = exp(-r*T);
terminalNoncarbonValue = vK*kappa-0.5*phiR*R.^2;
deterministicValue = runningValue-upgradeCostPV ...
    +discountTerminal*terminalNoncarbonValue;

maintenanceIntegral = trapz(tGrid,mLibrary,2);
maintenanceCarbonTotal = dt*sum(maintenanceEmissionRate,2);
upgradeCarbonTotal = sum(upgradeCarbonStep,2);
totalUpgradeIntensity = sum(omegaLibrary,2);
effectiveUpgradeCount = sum(omegaLibrary>0,2);

weightedUpgradeTime = nan(nPolicy,1);
for ic = 1:nPolicy
    if totalUpgradeIntensity(ic)>0
        eventTimeNow = eventStepLibrary(ic,:)*dt;
        weightedUpgradeTime(ic) = ...
            sum(eventTimeNow.*omegaLibrary(ic,:)) ...
            /totalUpgradeIntensity(ic);
    end
end

%% 5. Common random numbers for optimization and evaluation

% The same standard-normal shocks are rescaled by sigmaChi. This reduces
% Monte Carlo noise in comparisons across volatility levels.
% Antithetic shocks reduce finite-sample noise without changing the model.
halfOpt = floor(Mopt/2);
ZoptHalf = randn(halfOpt,N);
Zopt = [ZoptHalf;-ZoptHalf;zeros(Mopt-2*halfOpt,N)];

halfEval = floor(Meval/2);
ZevalHalf = randn(halfEval,N);
Zeval = [ZevalHalf;-ZevalHalf;zeros(Meval-2*halfEval,N)];

selectedPolicy = zeros(nSigma,1);
OmegaStar = zeros(nSigma,1);
upgradeCountStar = zeros(nSigma,1);
weightedTimeStar = nan(nSigma,1);
maintenanceStar = zeros(nSigma,1);

expectedObjective = zeros(nSigma,1);
objectiveCI95 = zeros(nSigma,1);

Pdynamic = zeros(nSigma,1);
Pterminal = zeros(nSigma,1);
Ptemporary = zeros(nSigma,1);
deficitSeverity = zeros(nSigma,1);

PcarbonRecovery = nan(nSigma,1);
expectedNetRecovery = nan(nSigma,1);

meanCarbonIntensity = zeros(nSigma,1);
sdCarbonIntensity = zeros(nSigma,1);

selectedMaintenancePath = zeros(nSigma,N+1);
selectedEventStep = zeros(nSigma,Lmax);
selectedOmega = zeros(nSigma,Lmax);

%% 6. Scenario-average optimization and independent evaluation

fprintf('\nStochastic carbon-intensity extension\n');
fprintf('Policies=%d, Mopt=%d, Meval=%d, N=%d\n\n', ...
    nPolicy,Mopt,Meval,N);

for is = 1:nSigma
    sigmaChi = sigmaChiGrid(is);

    % ---------- Optimization scenarios ----------
    chiOpt = chiBar*ones(Mopt,N+1);
    for n = 1:N
        chiOpt(:,n+1) = chiOpt(:,n) ...
            +muChi*(chiBar-chiOpt(:,n))*dt ...
            +sigmaChi*sqrt(dt)*Zopt(:,n);
        chiOpt(:,n+1) = min(max(chiOpt(:,n+1),chiLower),chiUpper);
    end
    chiOptMid = 0.5*(chiOpt(:,1:end-1)+chiOpt(:,2:end));

    % Scenario-specific terminal balance for every policy.
    operatingCarbonOpt = dt*(chiOptMid*energyCoefficient');
    terminalAOpt = A0+Q*T-operatingCarbonOpt ...
        -maintenanceCarbonTotal'-upgradeCarbonTotal';

    terminalSurplusOpt = max(terminalAOpt,0);
    terminalDeficitOpt = max(-terminalAOpt,0);
    terminalCarbonValueOpt = vA*terminalSurplusOpt ...
        -0.5*phiA*terminalDeficitOpt.^2;

    objectiveByPolicy = deterministicValue' ...
        +discountTerminal*mean(terminalCarbonValueOpt,1);

    [~,bestIndex] = max(objectiveByPolicy);
    selectedPolicy(is) = bestIndex;

    OmegaStar(is) = totalUpgradeIntensity(bestIndex);
    upgradeCountStar(is) = effectiveUpgradeCount(bestIndex);
    weightedTimeStar(is) = weightedUpgradeTime(bestIndex);
    maintenanceStar(is) = maintenanceIntegral(bestIndex);

    selectedMaintenancePath(is,:) = mLibrary(bestIndex,:);
    selectedEventStep(is,:) = eventStepLibrary(bestIndex,:);
    selectedOmega(is,:) = omegaLibrary(bestIndex,:);

    % ---------- Independent out-of-sample evaluation ----------
    chiEval = chiBar*ones(Meval,N+1);
    for n = 1:N
        chiEval(:,n+1) = chiEval(:,n) ...
            +muChi*(chiBar-chiEval(:,n))*dt ...
            +sigmaChi*sqrt(dt)*Zeval(:,n);
        chiEval(:,n+1) = min(max(chiEval(:,n+1),chiLower),chiUpper);
    end
    chiEvalMid = 0.5*(chiEval(:,1:end-1)+chiEval(:,2:end));

    % Realized moments are reported to verify that sigmaChi changes the
    % dispersion of chi without mechanically changing its sample mean.
    meanCarbonIntensity(is) = mean(chiEval(:));
    sdCarbonIntensity(is) = std(chiEval(:),0);

    Aeval = A0*ones(Meval,1);
    minimumAeval = Aeval;

    for n = 1:N
        Aeval = Aeval+dt*(Q ...
            -chiEvalMid(:,n)*energyCoefficient(bestIndex,n) ...
            -maintenanceEmissionRate(bestIndex,n)) ...
            -upgradeCarbonStep(bestIndex,n);
        minimumAeval = min(minimumAeval,Aeval);
    end

    terminalSurplusEval = max(Aeval,0);
    terminalDeficitEval = max(-Aeval,0);
    terminalCarbonValueEval = vA*terminalSurplusEval ...
        -0.5*phiA*terminalDeficitEval.^2;

    objectivePath = deterministicValue(bestIndex) ...
        +discountTerminal*terminalCarbonValueEval;

    expectedObjective(is) = mean(objectivePath);
    objectiveCI95(is) = 1.96*std(objectivePath,0)/sqrt(Meval);

    dynamicEvent = minimumAeval<0;
    terminalEvent = Aeval<0;
    temporaryEvent = dynamicEvent & ~terminalEvent;

    Pdynamic(is) = mean(dynamicEvent);
    Pterminal(is) = mean(terminalEvent);
    Ptemporary(is) = mean(temporaryEvent);
    deficitSeverity(is) = mean(max(-minimumAeval,0));

    % ---------- Portfolio carbon-recovery probability ----------
    % Counterfactual branch: same maintenance path, but remove all upgrades.
    if totalUpgradeIntensity(bestIndex)>0
        mSelected = mLibrary(bestIndex,:);
        kappaReference = kappa0;
        referenceEnergyCoefficient = zeros(1,N);

        for n = 1:N
            mMid = 0.5*(mSelected(n)+mSelected(n+1));
            kappaOld = kappaReference;
            kappaContinuous = kappaOld+dt*( ...
                etaK*mMid*(kappaBar-kappaOld)-deltaK*kappaOld);
            kappaContinuous = min(max(kappaContinuous,1e-8),kappaBar);
            kappaAverage = max(0.5*(kappaOld+kappaContinuous),1e-8);
            referenceEnergyCoefficient(n) = y/kappaAverage;
            kappaReference = kappaContinuous;
        end

        savingCoefficient = referenceEnergyCoefficient ...
            -energyCoefficient(bestIndex,:);
        operatingSavingEval = dt*sum(bsxfun(@times, ...
            chiEvalMid,savingCoefficient),2);
        netRecoveryEval = operatingSavingEval ...
            -upgradeCarbonTotal(bestIndex);

        PcarbonRecovery(is) = mean(netRecoveryEval>=0);
        expectedNetRecovery(is) = mean(netRecoveryEval);
    end

    fprintf(['sigmaChi=%4.2f | meanChi=%6.4f | sdChi=%6.4f | ' ...
        'Omega=%5.3f | M=%5.3f | E[J]=%8.4f | ' ...
        'Pdyn=%5.3f | PT=%5.3f | Ptemp=%5.3f\n'], ...
        sigmaChi,meanCarbonIntensity(is),sdCarbonIntensity(is), ...
        OmegaStar(is),maintenanceStar(is), ...
        expectedObjective(is),Pdynamic(is),Pterminal(is), ...
        Ptemporary(is));
end

%% 7. Results table and internal consistency checks

if any(Pterminal>Pdynamic+1e-12)
    error('Internal check failed: terminal deficit probability exceeds dynamic deficit probability.');
end

probabilityIdentityError = max(abs(Pdynamic-Pterminal-Ptemporary));
if probabilityIdentityError>1/Meval+1e-12
    warning('Pdynamic, Pterminal, and Ptemporary fail the sample identity check.');
end

% Mean-preservation diagnostic. With antithetic shocks and symmetric
% bounds, the error should be at machine precision.
meanPreservationError = max(abs(meanCarbonIntensity-chiBar));
if meanPreservationError>1e-10
    warning(['The realized mean carbon intensity differs from chiBar by ' ...
        '%.3e. Check the shock pairing and symmetric bounds.'], ...
        meanPreservationError);
end

% A boundary solution is not an error, but it must not be interpreted as
% an interior comparative-static result without an Lmax robustness check.
upgradeCapHit = OmegaStar>=Lmax-1e-8;
if any(upgradeCapHit)
    warning(['The total-upgrade upper bound is reached at one or more ' ...
        'volatility levels. Report an Lmax robustness check in the appendix.']);
end

% At zero volatility all scenarios coincide, so the confidence interval
% should be numerically zero. This is an internal benchmark check.
zeroVolatilityIndex = find(abs(sigmaChiGrid)<1e-14,1,'first');
if ~isempty(zeroVolatilityIndex) ...
        && objectiveCI95(zeroVolatilityIndex)>1e-10
    warning('The zero-volatility objective should have a zero sampling CI.');
end

ResultTable = table(sigmaChiGrid,OmegaStar,upgradeCountStar, ...
    weightedTimeStar,maintenanceStar,expectedObjective,objectiveCI95, ...
    Pdynamic,Pterminal,Ptemporary,deficitSeverity,PcarbonRecovery, ...
    expectedNetRecovery,meanCarbonIntensity,sdCarbonIntensity, ...
    upgradeCapHit,selectedPolicy, ...
    'VariableNames',{'SigmaChi','TotalUpgradeIntensity', ...
    'EffectiveUpgradeCount','WeightedUpgradeTime', ...
    'MaintenanceIntegral','ExpectedObjective','ObjectiveCI95', ...
    'P_DynamicDeficit','P_TerminalDeficit','P_TemporaryDeficit', ...
    'ExpectedMaximumDeficit','P_PortfolioCarbonRecovery', ...
    'ExpectedNetCarbonRecovery','MeanCarbonIntensity', ...
    'SDCarbonIntensity','UpgradeCapHit','SelectedPolicyIndex'});

writetable(ResultTable,resultCsv);

save(resultMat,'ResultTable','sigmaChiGrid','selectedPolicy', ...
    'selectedMaintenancePath','selectedEventStep','selectedOmega', ...
    'mLibrary','eventStepLibrary','omegaLibrary', ...
    'Mopt','Meval','N','dt','muChi','chiBar','chiLower','chiUpper', ...
    'meanPreservationError','upgradeCapHit');

%% 8. Publication-style figure

set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',11);

figure('Name','Figure 7: Stochastic carbon intensity extension', ...
    'Color','w','Position',[90 70 960 760]);

tl = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% (a) Total upgrade intensity
ax1 = nexttile;
plot(ax1,sigmaChiGrid,OmegaStar,'-o','Color','k', ...
    'LineWidth',1.7,'MarkerSize',5.5,'MarkerFaceColor',[0.25 0.25 0.25]);
xlabel(ax1,'Carbon-intensity volatility, $\sigma_{\chi}$', ...
    'Interpreter','latex');
ylabel(ax1,'Total upgrade intensity, $\Omega^*$', ...
    'Interpreter','latex');
title(ax1,'(a) Major upgrading');
grid(ax1,'on'); box(ax1,'on');

% (b) Cumulative maintenance effort
ax2 = nexttile;
plot(ax2,sigmaChiGrid,maintenanceStar,'-s','Color',[0.25 0.25 0.25], ...
    'LineWidth',1.7,'MarkerSize',5.5,'MarkerFaceColor','w');
xlabel(ax2,'Carbon-intensity volatility, $\sigma_{\chi}$', ...
    'Interpreter','latex');
ylabel(ax2,'Maintenance effort, $\int_0^T m^*(t)\,dt$', ...
    'Interpreter','latex');
title(ax2,'(b) Continuous maintenance');
grid(ax2,'on'); box(ax2,'on');

% (c) Out-of-sample economic performance
ax3 = nexttile;
errorbar(ax3,sigmaChiGrid,expectedObjective,objectiveCI95,'-o', ...
    'Color','k','LineWidth',1.6,'MarkerSize',5.5, ...
    'MarkerFaceColor',[0.25 0.25 0.25],'CapSize',7);
xlabel(ax3,'Carbon-intensity volatility, $\sigma_{\chi}$', ...
    'Interpreter','latex');
ylabel(ax3,'Out-of-sample expected objective');
title(ax3,'(c) Economic performance');
grid(ax3,'on'); box(ax3,'on');

% (d) Carbon-account downside risk
ax4 = nexttile;
h1 = plot(ax4,sigmaChiGrid,Pdynamic,'-o','Color','k', ...
    'LineWidth',1.7,'MarkerSize',5.5,'MarkerFaceColor',[0.25 0.25 0.25]);
hold(ax4,'on');
h2 = plot(ax4,sigmaChiGrid,Pterminal,'--s','Color',[0.35 0.35 0.35], ...
    'LineWidth',1.6,'MarkerSize',5.5,'MarkerFaceColor','w');
h3 = plot(ax4,sigmaChiGrid,Ptemporary,':^','Color',[0.55 0.55 0.55], ...
    'LineWidth',1.8,'MarkerSize',6.0,'MarkerFaceColor','w');
hold(ax4,'off');
xlabel(ax4,'Carbon-intensity volatility, $\sigma_{\chi}$', ...
    'Interpreter','latex');
ylabel(ax4,'Deficit probability');
title(ax4,'(d) Carbon-account downside risk');
ylim(ax4,[0 1]);
legend(ax4,[h1 h2 h3], ...
    {'$P_{\rm dyn}$','$P_T$', ...
    '$P_{\rm temp}=P_{\rm dyn}-P_T$'}, ...
    'Interpreter','latex','Location','northwest','Box','off');
grid(ax4,'on'); box(ax4,'on');

allAxes = [ax1 ax2 ax3 ax4];
for ia = 1:numel(allAxes)
    set(allAxes(ia),'LineWidth',1.0,'FontName','Times New Roman', ...
        'FontSize',11,'XLim',[min(sigmaChiGrid) max(sigmaChiGrid)]);
end


if exportFigures
    exportgraphics(gcf,figurePdf,'ContentType','vector');
    exportgraphics(gcf,figurePng,'Resolution',600);
end

%% 9. Console summary

fprintf('\nMain results saved to %s and %s.\n',resultCsv,resultMat);
fprintf('Probability identity max error = %.3e\n',probabilityIdentityError);
fprintf('Mean-intensity preservation max error = %.3e\n', ...
    meanPreservationError);
fprintf('Elapsed time = %.2f seconds.\n',toc);

fprintf(['\nInterpretation boundary:\n' ...
    '1) The policy is open-loop and common to all scenarios.\n' ...
    '2) The policy library approximates, rather than replaces, the ' ...
    'continuous-time optimality conditions.\n' ...
    '3) P_dynamic >= P_terminal by event inclusion. The informative ' ...
    'additional quantity is P_temporary.\n' ...
    '4) P_PortfolioCarbonRecovery compares the selected upgrade ' ...
    'portfolio with a no-upgrade branch under the same maintenance path.\n' ...
    '5) Symmetric winsorization preserves mean chi and prevents the ' ...
    'one-sided positivity correction from being mistaken for volatility.\n']);