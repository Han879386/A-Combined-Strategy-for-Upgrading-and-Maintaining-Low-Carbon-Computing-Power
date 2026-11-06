function Fig3_Fig4_StrictDynamicMaps_K0_Editable
%% ========================================================================
% Figure 3 and Figure 4 under the revised impulse-control formulation
%
% Figure 3:
%   Optimal maintenance regime in the continuation region.
%
% Figure 4:
%   Strict immediate-upgrade regime based on value matching:
%
% Output revision:
%   - No title is drawn inside either figure.
%   - MATLAB menus, toolbars and plot-edit interactions remain enabled.
%   - Editable .fig files are saved with the PDF and PNG files.
%       Delta_u(t,x) = M V(t,x) - C(t,x),
%   where C is the value of continuing without upgrading now and M V is the
%   best value of an effective upgrade now.
%
% Revisions relative to the old local-threshold code:
%   1) No representative shadow prices are imposed.
%   2) The fixed upgrade cost K0 enters the intervention comparison.
%   3) An active upgrade satisfies omega in [omegaMin,1].
%   4) The minimum upgrade interval is represented by a cooldown state.
%   5) No artificial upgrade-count cap is imposed.
%
% Preferred use:
%   Save the baseline dual-policy result from the path simulation:
%
%       save(fullfile(p.outputDir,'BaselineDualResult.mat'), ...
%            'res_dual','p','tGrid');
%
%   The code will then select a representative pre-upgrade state.
%% ========================================================================

clear; clc; close all;

p = baseline_parameters();
cfg = figure_settings();
slice = resolve_slice(cfg,p);

fprintf('\n================ Strict threshold-map slice ================\n');
fprintf('Requested t slice     = %.6f\n',slice.tRequested);
fprintf('Requested A slice     = %.6f\n',slice.ARequested);
if isfinite(slice.markerR)
    fprintf('Baseline marker R^-   = %.6f\n',slice.markerR);
    fprintf('Baseline marker G^-   = %.6f\n',slice.markerG);
end

switch lower(cfg.runMode)
    case 'fast'
        num.Nt = 50;
        num.nR = 19;
        num.nK = 19;
        num.nA = 21;
        num.nM = 9;
        num.nW = 17;
    case 'medium'
        num.Nt = 100;
        num.nR = 23;
        num.nK = 23;
        num.nA = 25;
        num.nM = 13;
        num.nW = 25;
    case 'final'
        num.Nt = 100;
        num.nR = 31;
        num.nK = 31;
        num.nA = 31;
        num.nM = 21;
        num.nW = 41;
    otherwise
        error('cfg.runMode must be fast, medium, or final.');
end

num.Rmax = 1.20;
num.Kmin = 0.25;
num.Amin = -1.50;
num.Amax = 3.20;
num.valueTol  = 1e-6;
num.actionTol = 2e-2;

dt = p.T/num.Nt;
disc = exp(-p.r*dt);
tGridDP = linspace(0,p.T,num.Nt+1);

gapRatio = p.minGap/dt;
gapSteps = round(gapRatio);
if abs(gapRatio-gapSteps)>1e-10
    warning(['minGap/dt is not an integer. The implemented interval is ' ...
        'rounded to %.6f.'],gapSteps*dt);
end

sliceIndex = find(tGridDP(1:end-1)>=slice.tRequested-1e-12,1,'first');
if isempty(sliceIndex)
    error('The requested time slice is outside the decision horizon.');
end
actualTime = tGridDP(sliceIndex);

Rgrid = make_grid_with_anchor(0,num.Rmax,num.nR, ...
    finite_anchor(slice.markerR,p.R0));
Kanchor = p.kappaBar - finite_anchor(slice.markerG, ...
    p.kappaBar-p.kappa0);
Kgrid = make_grid_with_anchor(num.Kmin,p.kappaBar,num.nK,Kanchor);
Agrid = make_grid_with_anchor(num.Amin,num.Amax,num.nA,slice.ARequested);

[R3,K3,A3] = ndgrid(Rgrid,Kgrid,Agrid);

mSet = linspace(0,p.mbar,num.nM);
wSet = linspace(p.omegaMin,1,num.nW);

[~,aIndex] = min(abs(Agrid-slice.ARequested));
actualA = Agrid(aIndex);

maintTrans = cell(num.nM,1);
for im=1:num.nM
    m=mSet(im);

    Rn = R3 + dt.*(p.beta*p.y-(p.deltaR+p.etaR*m).*R3);
    Kn = K3 + dt.*(p.etaK*m.*(p.kappaBar-K3)-p.deltaK.*K3);
    An = A3 + dt.*(p.Q-p.chi*p.y./max(K3,1e-8) ...
        -p.em*m-0.5*p.em2*m.^2-p.deltaA.*A3);

    [Rn,Kn,An]=clip_states(Rn,Kn,An,Rgrid,Kgrid,Agrid);

    profit = p.price*p.y ...
        -p.pE*p.y./max(K3,1e-8) ...
        -0.5*p.ay.*(p.y./max(K3,1e-8)).^2 ...
        -p.cm*m-0.5*p.am*m.^2 ...
        -0.5*p.aR.*R3.^2;

    maintTrans{im}=struct( ...
        'R',Rn,'K',Kn,'A',An,'reward',dt.*profit,'m',m);
end

upgradeTrans = cell(num.nW,1);
for iw=1:num.nW
    w=wSet(iw);

    Rp=(1-p.gammaR*w).*R3;
    Kp=K3+p.gammaK*w.*(p.kappaBar-K3);
    Ap=A3-p.xi1*w-0.5*p.xi2*w.^2;

    [Rp,Kp,Ap]=clip_states(Rp,Kp,Ap,Rgrid,Kgrid,Agrid);

    upgradeTrans{iw}=struct( ...
        'R',Rp,'K',Kp,'A',Ap, ...
        'cost',p.K0+p.K1*w+0.5*p.K2*w.^2, ...
        'w',w);
end

terminalV = p.vA.*max(A3,0) ...
    -0.5*p.phiA.*max(-A3,0).^2 ...
    +p.vK.*K3 ...
    -0.5*p.phiR.*R3.^2;

Vnext = repmat(terminalV,[1,1,1,gapSteps+1]);

sliceSaved=false;
mSlice=[];
deltaSlice=[];
omegaSlice=[];
upgradeClassSlice=[];

fprintf('\nBackward dynamic program started (%s mode).\n',cfg.runMode);
fprintf('Grid: R=%d, kappa=%d, A=%d, time=%d, cooldown=%d\n', ...
    num.nR,num.nK,num.nA,num.Nt,gapSteps+1);

for n=num.Nt:-1:1
    timeNow=tGridDP(n);

    calendarEligible = ...
        timeNow>=p.sigmaMin-1e-12 && ...
        timeNow<=p.sigmaMax+1e-12;

    Ccell=cell(gapSteps+1,1);
    MpolCell=cell(gapSteps+1,1);

    for c=0:gapSteps
        cNext=max(c-1,0);

        Fnext=griddedInterpolant( ...
            {Rgrid,Kgrid,Agrid},Vnext(:,:,:,cNext+1), ...
            'linear','nearest');

        Cbest=-inf(size(R3));
        mbest=zeros(size(R3));

        for im=1:num.nM
            tr=maintTrans{im};
            candidate=tr.reward+disc.*Fnext(tr.R,tr.K,tr.A);

            improve=candidate>Cbest;
            Cbest(improve)=candidate(improve);
            mbest(improve)=tr.m;
        end

        Ccell{c+1}=Cbest;
        MpolCell{c+1}=mbest;
    end

    Vcurr=zeros(size(Vnext));

    for c=0:gapSteps
        Cbest=Ccell{c+1};
        Mbest=-inf(size(R3));
        wbest=zeros(size(R3));

        if c==0 && calendarEligible
            Fpost=griddedInterpolant( ...
                {Rgrid,Kgrid,Agrid},Ccell{gapSteps+1}, ...
                'linear','nearest');

            for iw=1:num.nW
                trU=upgradeTrans{iw};
                candidate=Fpost(trU.R,trU.K,trU.A)-trU.cost;

                improve=candidate>Mbest;
                Mbest(improve)=candidate(improve);
                wbest(improve)=trU.w;
            end
        end

        Vcurr(:,:,:,c+1)=max(Cbest,Mbest);

        if n==sliceIndex && c==0
            deltaFull=Mbest-Cbest;
            mFull=MpolCell{1};

            mSlice=squeeze(mFull(:,:,aIndex));
            deltaSlice=squeeze(deltaFull(:,:,aIndex));
            omegaSlice=squeeze(wbest(:,:,aIndex));

            upgradeClassSlice=zeros(size(deltaSlice));
            active=deltaSlice>num.valueTol;
            partial=active & omegaSlice<1-num.actionTol;
            full=active & omegaSlice>=1-num.actionTol;

            upgradeClassSlice(partial)=1;
            upgradeClassSlice(full)=2;
            sliceSaved=true;
        end
    end

    Vnext=Vcurr;

    if mod(num.Nt-n+1,10)==0 || n==1
        fprintf('  completed %d/%d backward steps\n', ...
            num.Nt-n+1,num.Nt);
    end
end

if ~sliceSaved
    error('The requested threshold slice was not saved.');
end

maintenanceClass=ones(size(mSlice));
maintenanceClass(mSlice<=num.actionTol)=0;
maintenanceClass(mSlice>=p.mbar-num.actionTol)=2;

upgradeActive=upgradeClassSlice>0;
maintenanceMasked=double(maintenanceClass);
maintenanceMasked(upgradeActive)=NaN;

Graw=p.kappaBar-Kgrid;
[Gplot,gOrder]=sort(Graw,'ascend');

maintenancePlot=maintenanceMasked(:,gOrder)';
upgradePlot=upgradeClassSlice(:,gOrder)';
deltaPlot=deltaSlice(:,gOrder)';
omegaPlot=omegaSlice(:,gOrder)';
activePlot=upgradeActive(:,gOrder)';

set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',9.5);

cmap=[0.88 0.91 0.94;0.64 0.70 0.75;0.30 0.34 0.38];

fig3=figure('Name','','NumberTitle','off', ...
    'Color','w','Units','centimeters','Position',[3,3,14.8,11.2], ...
    'Renderer','painters','MenuBar','figure','ToolBar','figure');
ax3=axes(fig3);

h3=imagesc(ax3,Rgrid,Gplot,maintenancePlot);
axis(ax3,'xy');
set(h3,'AlphaData',~isnan(maintenancePlot));
ax3.Color=[0.95 0.95 0.95];
colormap(ax3,cmap);
caxis(ax3,[-0.5 2.5]);
hold(ax3,'on');

contour(ax3,Rgrid,Gplot,maintenancePlot,[0.5 1.5], ...
    'k','LineWidth',1.1);

if any(activePlot(:)) && any(~activePlot(:))
    contour(ax3,Rgrid,Gplot,deltaPlot,[0 0], ...
        'k--','LineWidth',1.4);
end

plot_reference_marker(ax3,slice);
hold(ax3,'off');

xlabel(ax3,'Privacy risk, $R$','Interpreter','latex');
ylabel(ax3,'Efficiency gap, $G=\bar{\kappa}-\kappa$', ...
    'Interpreter','latex');

cb3=colorbar(ax3);
cb3.Ticks=[0 1 2];
cb3.TickLabels={'No maintenance','Interior','Maximum'};
style_axis(ax3);

if cfg.exportFigures
    savefig(fig3,'Figure3_StrictMaintenanceMap_K0_Editable.fig');
    exportgraphics(fig3,'Figure3_StrictMaintenanceMap_K0.pdf', ...
        'ContentType','vector');
    exportgraphics(fig3,'Figure3_StrictMaintenanceMap_K0.png', ...
        'Resolution',600);
end

fig4=figure('Name','','NumberTitle','off', ...
    'Color','w','Units','centimeters','Position',[4,4,14.8,11.2], ...
    'Renderer','painters','MenuBar','figure','ToolBar','figure');
ax4=axes(fig4);

imagesc(ax4,Rgrid,Gplot,upgradePlot);
axis(ax4,'xy');
colormap(ax4,cmap);
caxis(ax4,[-0.5 2.5]);
hold(ax4,'on');

finiteDelta=deltaPlot(isfinite(deltaPlot));
if ~isempty(finiteDelta) && min(finiteDelta)<0 && max(finiteDelta)>0
    contour(ax4,Rgrid,Gplot,deltaPlot,[0 0], ...
        'k--','LineWidth',1.5);
end

if any(upgradePlot(:)==1) && any(upgradePlot(:)==2)
    contour(ax4,Rgrid,Gplot,omegaPlot, ...
        [1-num.actionTol 1-num.actionTol], ...
        'k','LineWidth',1.1);
end

plot_reference_marker(ax4,slice);
hold(ax4,'off');

xlabel(ax4,'Pre-upgrade privacy risk, $R(\sigma_l^-)$', ...
    'Interpreter','latex');
ylabel(ax4,'Pre-upgrade efficiency gap, $G(\sigma_l^-)$', ...
    'Interpreter','latex');

cb4=colorbar(ax4);
cb4.Ticks=[0 1 2];
cb4.TickLabels={'No upgrade','Partial upgrade','Full upgrade'};
style_axis(ax4);

if cfg.exportFigures
    savefig(fig4,'Figure4_StrictUpgradeMap_K0_Editable.fig');
    exportgraphics(fig4,'Figure4_StrictUpgradeMap_K0.pdf', ...
        'ContentType','vector');
    exportgraphics(fig4,'Figure4_StrictUpgradeMap_K0.png', ...
        'Resolution',600);
end

fprintf('\n================ Threshold-map diagnostics ================\n');
fprintf('Displayed time                         : %.6f\n',actualTime);
fprintf('Displayed carbon balance               : %.6f\n',actualA);
fprintf('Implemented minimum interval           : %.6f\n',gapSteps*dt);
fprintf('Maintenance no/interior/maximum        : %d / %d / %d\n', ...
    nnz(maintenanceClass==0 & ~upgradeActive), ...
    nnz(maintenanceClass==1 & ~upgradeActive), ...
    nnz(maintenanceClass==2 & ~upgradeActive));
fprintf('Immediate-upgrade states masked in Fig3: %d\n', ...
    nnz(upgradeActive));
fprintf('Upgrade no/partial/full                : %d / %d / %d\n', ...
    nnz(upgradeClassSlice==0), ...
    nnz(upgradeClassSlice==1), ...
    nnz(upgradeClassSlice==2));

if ~isempty(finiteDelta)
    fprintf('Delta_u range                          : [%.8f, %.8f]\n', ...
        min(finiteDelta),max(finiteDelta));
end

% Keep the axes interactive and open both figures in plot-edit mode.
try
    enableDefaultInteractivity(ax3);
    enableDefaultInteractivity(ax4);
    plotedit(fig3,'on');
    plotedit(fig4,'on');
catch
    warning(['Automatic plot-edit mode is unavailable in this MATLAB ' ...
        'version. Reopen the saved .fig file and execute plotedit on.']);
end

end

function cfg=figure_settings()
cfg.runMode='fast';       % fast | medium | final
cfg.exportFigures=true;
cfg.useBaselineResult=true;
cfg.baselineResultFile=fullfile('publication_figures', ...
    'BaselineDualResult.mat');
cfg.eventIndex=[];
cfg.manualTime=3.60;
cfg.manualA=1.50;
end

function p=baseline_parameters()
p.T=10;
p.r=0.05;

p.y=1.00;
p.price=3.00;
p.pE=0.60;
p.chi=0.35;
p.Q=0.50;
p.deltaA=0.02;

p.R0=0.40;
p.kappa0=0.55;
p.A0=2.20;
p.kappaBar=1.00;

p.beta=0.25;
p.deltaR=0.05;
p.etaR=0.45;

p.etaK=0.35;
p.deltaK=0.12;

p.cm=0.08;
p.am=0.45;
p.em=0.03;
p.em2=0.10;

p.aR=1.20;
p.ay=0.45;

p.K0=0.05;
p.K1=0.25;
p.K2=0.70;

p.xi1=0.14;
p.xi2=0.35;

p.gammaR=0.65;
p.gammaK=0.65;

p.vA=0.50;
p.phiA=8.00;
p.vK=0.80;
p.phiR=1.50;

p.mbar=1.00;

p.sigmaMin=0.50;
p.sigmaMax=9.50;
p.minGap=0.80;
p.omegaMin=0.03;
end

function slice=resolve_slice(cfg,p)
slice.tRequested=cfg.manualTime;
slice.ARequested=cfg.manualA;
slice.markerR=NaN;
slice.markerG=NaN;

if ~cfg.useBaselineResult || ~isfile(cfg.baselineResultFile)
    if cfg.useBaselineResult
        warning(['BaselineDualResult.mat was not found. The manual slice ' ...
            '(t=%.2f,A=%.2f) is used.'], ...
            cfg.manualTime,cfg.manualA);
    end
    return;
end

S=load(cfg.baselineResultFile);
if ~isfield(S,'res_dual')
    warning('The MAT file does not contain res_dual. Manual slice is used.');
    return;
end

res=S.res_dual;
L=numel(res.sigma);
if L==0
    warning('The baseline policy has no active upgrade. Manual slice is used.');
    return;
end

if isempty(cfg.eventIndex)
    lRef=ceil(L/2);
else
    lRef=min(max(round(cfg.eventIndex),1),L);
end

tRef=res.sigma(lRef);
idx=find(abs(res.pathTime-tRef)<1e-7,1,'first');
if isempty(idx)
    [~,idx]=min(abs(res.pathTime-tRef));
end

slice.tRequested=tRef;
slice.ARequested=res.Apath(idx);
slice.markerR=res.Rpath(idx);
slice.markerG=p.kappaBar-res.Kpath(idx);

fprintf('Baseline event %d of %d selected for the threshold slice.\n', ...
    lRef,L);
end

function anchor=finite_anchor(candidate,fallback)
if isfinite(candidate)
    anchor=candidate;
else
    anchor=fallback;
end
end

function grid=make_grid_with_anchor(lo,hi,n,anchor)
grid=linspace(lo,hi,n);
anchor=min(max(anchor,lo),hi);
[~,idx]=min(abs(grid-anchor));
grid(idx)=anchor;
grid=sort(grid);
if any(diff(grid)<=0)
    grid=unique([linspace(lo,hi,n),anchor]);
end
end

function [R,K,A]=clip_states(R,K,A,Rgrid,Kgrid,Agrid)
R=min(max(R,Rgrid(1)),Rgrid(end));
K=min(max(K,Kgrid(1)),Kgrid(end));
A=min(max(A,Agrid(1)),Agrid(end));
end

function plot_reference_marker(ax,slice)
if isfinite(slice.markerR) && isfinite(slice.markerG)
    scatter(ax,slice.markerR,slice.markerG,58,'p', ...
        'MarkerFaceColor','w','MarkerEdgeColor','k', ...
        'LineWidth',1.0,'DisplayName','Baseline pre-upgrade state');
end
end

function style_axis(ax)
grid(ax,'off');
box(ax,'on');
ax.Layer='top';
ax.LineWidth=0.85;
ax.FontName='Times New Roman';
ax.FontSize=9.5;
ax.TickDir='out';
ax.XMinorTick='on';
ax.YMinorTick='on';
end

