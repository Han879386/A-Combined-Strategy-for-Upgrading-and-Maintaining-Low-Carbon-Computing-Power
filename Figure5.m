%% ========================================================================
%  Figure 5: Carbon payback of a major upgrade
%
%  TRUE MATLAB SCRIPT: no user-defined or local functions are included.
%  Save this file as Figure5_AccountCarbonPayback_Proposition3.m and click Run.
%
%  Purpose
%    (a) Compare account-adjusted net carbon recovery on CALENDAR TIME.
%    (b) Compare latest upgrade times that complete ACCOUNT carbon payback.
%    Both panels directly implement the revised Proposition 3.
%
%  Identification
%    For every candidate upgrade, the reference and upgrade branches start
%    from the same pre-upgrade state and use the same subsequent maintenance
%    path. Thus, their emission difference isolates the carbon effect of the
%    major upgrade. Maintenance emissions and allowance inflows cancel in
%    the branch difference.
%
%  Carbon-account model
%    dA/dt = Q - chi*y/kappa - em*m - 0.5*em2*m^2 - deltaA*A
%    E_U(omega) = xi1*omega + 0.5*xi2*omega^2
%
%  Definitions for an upgrade at time sigma, elapsed horizon u:
%    Physical net saving: integral_0^u Deltae(v) dv - E_U
%    Adjusted net saving: integral_0^u exp(deltaA*v)*Deltae(v) dv - E_U
%    Terminal account difference: exp(-deltaA*(T-sigma))*adjusted_net(T)
%
%  This is a fixed-continuation-policy COUNTERFACTUAL comparison, not
%  the globally optimal upgrade schedule and NOT a pathwise safety test.
%  Maintenance and allowance inflows cancel between both branches.
%% ========================================================================

clearvars;
clc;
close all;

%% =========================== Baseline parameters ========================

p.T = 10;
p.r = 0.05;

p.y = 1.00;
p.price = 3.00;
p.pE = 0.60;
p.chi = 0.35;
p.Q = 0.50;
p.deltaA = 0.02;  % carbon-account adjustment rate (baseline manuscript)

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

p.K0 = 0.08;
p.K1 = 0.25;
p.K2 = 0.70;

p.xi1 = 0.14;
p.xi2 = 0.35;

p.gammaR = 0.65;
p.gammaK = 0.70;

p.vA = 0.50;
p.phiA = 8.00;
p.vK = 0.80;
p.phiR = 1.50;

p.mbar = 1.00;

%% =========================== Analysis settings ==========================

% A fine fixed time step is sufficient because this calculation contains
% no optimization. All candidate times below are exact grid points.
dt = 0.01;
tGrid = (0:dt:p.T)';
N = numel(tGrid)-1;

% Local-envelope comparison: use the same moderate continuation
% maintenance effort in the reference and upgrade branches.
mBenchmark = 0.50*p.mbar;
mPath = mBenchmark*ones(size(tGrid));

% Candidate upgrade times and intensities.
tCandidate = (0.50:0.05:9.50)';
omegaSet = [0.40 0.60 0.80];

% Panel (a) uses the benchmark intensity and two timing cases.
omegaBenchmark = 0.60;
tEarlier = 3.00;
tLater = 8.00;

exportFigures = true;
outputFig = 'Figure5_AccountCarbonPayback_Proposition3.fig';
outputPdf = 'Figure5_AccountCarbonPayback_Proposition3.pdf';
outputPng = 'Figure5_AccountCarbonPayback_Proposition3.png';
outputMat = 'Figure5_AccountCarbonPayback_Proposition3.mat';

%% ======================= Pre-upgrade benchmark path ====================

% The pre-upgrade state follows the model without any previous impulse.
% Only kappa is required because both policies share other carbon inputs.
kappaBase = zeros(size(tGrid));
kappaBase(1) = p.kappa0;

for n = 1:N
    kOld = kappaBase(n);
    mNow = mPath(n);

    kappaBase(n+1) = kOld+dt*( ...
        p.etaK*mNow*(p.kappaBar-kOld)-p.deltaK*kOld);

    kappaBase(n+1) = min(max(kappaBase(n+1),1e-8),p.kappaBar);
end

%% ==================== Carbon-payback calculations ======================

nTime = numel(tCandidate);
nOmega = numel(omegaSet);

paybackTime = nan(nTime,nOmega);  % account-adjusted recovery time
recoveryRatio = zeros(nTime,nOmega); % ACCOUNT recovery ratio
physicalRecoveryRatio = zeros(nTime,nOmega);
remainingHorizon = p.T-tCandidate;
terminalCarbonAdvantage = zeros(nTime,nOmega); % true terminal A_U-A_0
terminalPhysicalNetSaving = zeros(nTime,nOmega);
terminalAccountNetRecovery = zeros(nTime,nOmega);

elapsedEarlier = [];
elapsedLater = [];
calendarEarlier = [];
calendarLater = [];
netCarbonEarlier = [];
netCarbonLater = [];
paybackEarlier = NaN;
paybackLater = NaN;

for iw = 1:nOmega
    omega = omegaSet(iw);
    upgradeCarbon = p.xi1*omega+0.5*p.xi2*omega^2;

    for it = 1:nTime
        sigma = tCandidate(it);
        eventIndex = round(sigma/dt)+1;
        branchTime = tGrid(eventIndex:end);
        elapsed = branchTime-sigma;
        nBranch = numel(branchTime);

        kappaReference = zeros(nBranch,1);
        kappaUpgrade = zeros(nBranch,1);

        kappaMinus = kappaBase(eventIndex);
        kappaReference(1) = kappaMinus;
        kappaUpgrade(1) = kappaMinus ...
            +p.gammaK*omega*(p.kappaBar-kappaMinus);

        for j = 1:nBranch-1
            globalIndex = eventIndex+j-1;
            mNow = mPath(globalIndex);

            kRefOld = kappaReference(j);
            kUpOld = kappaUpgrade(j);

            kappaReference(j+1) = kRefOld+dt*( ...
                p.etaK*mNow*(p.kappaBar-kRefOld)-p.deltaK*kRefOld);

            kappaUpgrade(j+1) = kUpOld+dt*( ...
                p.etaK*mNow*(p.kappaBar-kUpOld)-p.deltaK*kUpOld);

            kappaReference(j+1) = min(max( ...
                kappaReference(j+1),1e-8),p.kappaBar);
            kappaUpgrade(j+1) = min(max( ...
                kappaUpgrade(j+1),1e-8),p.kappaBar);
        end

        mBranch = mPath(eventIndex:end);

        emissionReference = p.chi*p.y./kappaReference ...
            +p.em*mBranch+0.5*p.em2*mBranch.^2;

        emissionUpgrade = p.chi*p.y./kappaUpgrade ...
            +p.em*mBranch+0.5*p.em2*mBranch.^2;

        avoidedEmissionRate = emissionReference-emissionUpgrade;
        physicalOperatingSaving = cumtrapz(elapsed,avoidedEmissionRate);
        accountOperatingSaving = cumtrapz(elapsed, ...
            exp(p.deltaA*elapsed).*avoidedEmissionRate);

        netPhysicalCarbon = physicalOperatingSaving-upgradeCarbon;
        netCumulativeCarbon = accountOperatingSaving-upgradeCarbon;

        physicalRecoveryRatio(it,iw) = ...
            physicalOperatingSaving(end)/upgradeCarbon;
        recoveryRatio(it,iw) = ...
            accountOperatingSaving(end)/upgradeCarbon;
        terminalPhysicalNetSaving(it,iw) = netPhysicalCarbon(end);
        terminalAccountNetRecovery(it,iw) = netCumulativeCarbon(end);
        terminalCarbonAdvantage(it,iw) = ...
            exp(-p.deltaA*elapsed(end))*netCumulativeCarbon(end);

        crossingIndex = find(netCumulativeCarbon>=0,1,'first');

        if ~isempty(crossingIndex)
            if crossingIndex == 1
                tauC = 0;
            else
                x1 = elapsed(crossingIndex-1);
                x2 = elapsed(crossingIndex);
                y1 = netCumulativeCarbon(crossingIndex-1);
                y2 = netCumulativeCarbon(crossingIndex);

                if abs(y2-y1) < 1e-14
                    tauC = x2;
                else
                    tauC = x1-y1*(x2-x1)/(y2-y1);
                end
            end

            paybackTime(it,iw) = tauC;
        end

        if abs(omega-omegaBenchmark)<1e-12 ...
                && abs(sigma-tEarlier)<1e-12
            elapsedEarlier = elapsed;
            calendarEarlier = branchTime;
            netCarbonEarlier = netCumulativeCarbon;
            paybackEarlier = paybackTime(it,iw);
        end

        if abs(omega-omegaBenchmark)<1e-12 ...
                && abs(sigma-tLater)<1e-12
            elapsedLater = elapsed;
            calendarLater = branchTime;
            netCarbonLater = netCumulativeCarbon;
            paybackLater = paybackTime(it,iw);
        end
    end
end

% Interpolate latest candidate time satisfying ACCOUNT CRR_U>=1.
% The base-case CRR is strictly decreasing in upgrade time. Check this
% assumption before interpreting the interpolated boundary as unique.
latestFeasibleTime = nan(1,nOmega);

for iw = 1:nOmega
    ratioNow = recoveryRatio(:,iw);

    if any(diff(ratioNow)>1e-8)
        warning('Account CRR is not monotone for omega = %.2f; inspect candidate windows.', ...
            omegaSet(iw));
    end

    lastAbove = find(ratioNow>=1,1,'last');
    if isempty(lastAbove)
        latestFeasibleTime(iw) = NaN;
    elseif lastAbove == nTime
        latestFeasibleTime(iw) = tCandidate(end);
    else

        x1 = tCandidate(lastAbove);
        x2 = tCandidate(lastAbove+1);
        y1 = ratioNow(lastAbove);
        y2 = ratioNow(lastAbove+1);

        if abs(y2-y1)<1e-14
            latestFeasibleTime(iw) = x1;
        else
            latestFeasibleTime(iw) = ...
                x1+(1-y1)*(x2-x1)/(y2-y1);
        end
    end
end

%% ================================ Plot ==================================

set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',10.5);

fig5 = figure('Name','Figure 5: Account-adjusted payback of major upgrading', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[3,3,18.4,8.4], ...
    'Renderer','painters', ...
    'MenuBar','figure', ...
    'ToolBar','figure');

tl = tiledlayout(fig5,1,2, ...
    'Padding','compact','TileSpacing','compact');

% ------------------ (a) Account-adjusted net carbon recovery ------------
ax1 = nexttile(tl,1);
hold(ax1,'on');

hEarly = plot(ax1,calendarEarlier,netCarbonEarlier,'k-', ...
    'LineWidth',1.8);
hLate = plot(ax1,calendarLater,netCarbonLater,'--', ...
    'Color',[0.45 0.45 0.45], ...
    'LineWidth',1.8);
yline(ax1,0,'k:','LineWidth',1.1,'HandleVisibility','off');

if isfinite(paybackEarlier)
    calendarPaybackEarlier = tEarlier+paybackEarlier;
    plot(ax1,calendarPaybackEarlier,0,'ko', ...
        'MarkerFaceColor','k','MarkerSize',5.5, ...
        'HandleVisibility','off');
    text(ax1,calendarPaybackEarlier,0, ...
        sprintf('  %.2f',calendarPaybackEarlier), ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'FontSize',8.5);
end

if isfinite(paybackLater)
    calendarPaybackLater = tLater+paybackLater;
    plot(ax1,calendarPaybackLater,0,'ks', ...
        'MarkerFaceColor','w','MarkerSize',5.5, ...
        'HandleVisibility','off');
    text(ax1,calendarPaybackLater,0, ...
        sprintf('%.2f  ',calendarPaybackLater), ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','top', ...
        'FontSize',8.5);
end

hold(ax1,'off');
xlabel(ax1,'Calendar time, $t$', ...
    'Interpreter','latex');
ylabel(ax1,'$B_U(u,\omega;\sigma)$ (account-adjusted)', ...
    'Interpreter','latex');
title(ax1,sprintf('(a) Account carbon payback, $\\omega=%.1f$',omegaBenchmark), ...
    'Interpreter','latex','FontWeight','normal');
xlim(ax1,[2.5 p.T]);

legend(ax1,[hEarly hLate], ...
    {sprintf('Earlier upgrade, $t=%.1f$',tEarlier), ...
     sprintf('Later upgrade, $t=%.1f$',tLater)}, ...
    'Interpreter','latex','Location','southeast','Box','off');

grid(ax1,'on');
box(ax1,'on');

% ------------------ (b) Latest upgrade time completing payback ---------
ax2 = nexttile(tl,2);
hold(ax2,'on');

hLatest = plot(ax2,omegaSet,latestFeasibleTime,'ko-', ...
    'LineWidth',1.7, ...
    'MarkerSize',6.0, ...
    'MarkerFaceColor',[0.30 0.30 0.30]);

for iw = 1:nOmega
    if isfinite(latestFeasibleTime(iw))
        text(ax2,omegaSet(iw),latestFeasibleTime(iw), ...
            sprintf('  %.2f',latestFeasibleTime(iw)), ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','bottom', ...
            'FontSize',9.0);
    end
end

hold(ax2,'off');

xlabel(ax2,'Upgrade intensity, $\omega$', ...
    'Interpreter','latex');
ylabel(ax2,'Latest account-payback time, $\sigma_C^{\max}$', ...
    'Interpreter','latex');
title(ax2,'(b) Latest upgrade time allowing payback', ...
    'FontWeight','normal');

set(ax2,'XTick',omegaSet);
xlim(ax2,[min(omegaSet)-0.08 max(omegaSet)+0.12]);

finiteLatest = latestFeasibleTime(isfinite(latestFeasibleTime));
if ~isempty(finiteLatest)
    yPadding = max(0.15,0.35*(max(finiteLatest)-min(finiteLatest)+0.05));
    ylim(ax2,[min(finiteLatest)-yPadding max(finiteLatest)+yPadding]);
end

legend(ax2,hLatest,{'$\mathrm{CRR}_U^{acc}(\sigma_C^{\max},\omega)=1$'}, ...
    'Interpreter','latex','Location','southwest','Box','off');

grid(ax2,'on');
box(ax2,'on');

set([ax1 ax2], ...
    'LineWidth',0.9, ...
    'FontName','Times New Roman', ...
    'FontSize',10.5, ...
    'TickDir','out', ...
    'XMinorTick','on', ...
    'YMinorTick','on');

%% =============================== Diagnostics ============================

fprintf('\n================ Carbon-payback analysis ================\n');
fprintf('Planning horizon T                    = %.2f\n',p.T);
fprintf('Carbon account adjustment deltaA      = %.4f\n',p.deltaA);
fprintf('Benchmark maintenance effort          = %.4f\n',mBenchmark);
fprintf('Benchmark upgrade intensity           = %.4f\n',omegaBenchmark);
fprintf('One-off benchmark upgrade emissions   = %.6f\n', ...
    p.xi1*omegaBenchmark+0.5*p.xi2*omegaBenchmark^2);

fprintf('\nIntensity   tau_C(t=%.1f)   AccCRR(t=%.1f)   tau_C(t=%.1f)   AccCRR(t=%.1f)   Latest payback sigma\n', ...
    tEarlier,tEarlier,tLater,tLater);

for iw = 1:nOmega
    [~,iEarly] = min(abs(tCandidate-tEarlier));
    [~,iLate] = min(abs(tCandidate-tLater));

    latestFeasible = latestFeasibleTime(iw);

    fprintf('  %.2f        ',omegaSet(iw));

    if isfinite(paybackTime(iEarly,iw))
        fprintf('%9.4f      ',paybackTime(iEarly,iw));
    else
        fprintf('  no payback      ');
    end

    fprintf('%9.4f      ',recoveryRatio(iEarly,iw));

    if isfinite(paybackTime(iLate,iw))
        fprintf('%9.4f      ',paybackTime(iLate,iw));
    else
        fprintf('  no payback      ');
    end

    fprintf('%9.4f      ',recoveryRatio(iLate,iw));

    if isfinite(latestFeasible)
        fprintf('%9.4f\n',latestFeasible);
    else
        fprintf('none\n');
    end
end

fprintf(['\nInterpretation: AccCRR_U>=1 means the account-adjusted ' ...
    'future operating savings recover the initial upgrade emission ' ...
    'in the terminal carbon-account comparison (Proposition 3). ' ...
    'It is NOT a pathwise carbon-account safety check.\n']);
[~,iEarlyDiag] = min(abs(tCandidate-tEarlier));
for iw=1:nOmega
    fprintf(['omega %.2f: physical CRR at t=%.1f is %.4f; ' ...
        'account CRR is %.4f\n'], ...
        omegaSet(iw),tEarlier, ...
        physicalRecoveryRatio(iEarlyDiag,iw),recoveryRatio(iEarlyDiag,iw));
end

if isfinite(paybackEarlier)
    fprintf('Benchmark earlier-upgrade payback time = %.6f\n', ...
        paybackEarlier);
else
    fprintf('The benchmark earlier upgrade does not pay back by T.\n');
end

if isfinite(paybackLater)
    fprintf('Benchmark later-upgrade payback time   = %.6f\n', ...
        paybackLater);
else
    fprintf('The benchmark later upgrade does not pay back by T.\n');
end

%% ================================ Export ================================

if exportFigures
    savefig(fig5,outputFig);
    exportgraphics(fig5,outputPdf,'ContentType','vector');
    exportgraphics(fig5,outputPng,'Resolution',600);

    save(outputMat,'p','dt','tGrid','mBenchmark','mPath', ...
        'tCandidate','omegaSet','omegaBenchmark','tEarlier','tLater', ...
        'paybackTime','recoveryRatio','physicalRecoveryRatio','remainingHorizon', ...
        'latestFeasibleTime', ...
        'terminalCarbonAdvantage','terminalPhysicalNetSaving', ...
        'terminalAccountNetRecovery','elapsedEarlier','elapsedLater', ...
        'calendarEarlier','calendarLater', ...
        'netCarbonEarlier','netCarbonLater','paybackEarlier', ...
        'paybackLater','kappaBase');

    fprintf('\nSaved files:\n');
    fprintf('  %s\n',outputFig);
    fprintf('  %s\n',outputPdf);
    fprintf('  %s\n',outputPng);
    fprintf('  %s\n',outputMat);
end

try
    enableDefaultInteractivity(ax1);
    enableDefaultInteractivity(ax2);
catch
    % Older MATLAB versions may not provide this function.
end
