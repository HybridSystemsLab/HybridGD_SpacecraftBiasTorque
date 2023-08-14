%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
%--------------------------------------------------------------------------
% Project: Simulation of a hybrid system
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

clearvars
close all
clc

%% Parameters
Js = 5000;
Jw = 0.1;
delta = 9.5;
Mt = -10;
theta = 0.005;
u = 0;

Kp = 10;
Kd = 1200;
Ki = 0.00045;

Omega_max = convangvel(10000,'rpm','rad/s');
taus = 10;

% estimator
gammac = 0.0012;
lambdac = 0.001; % 0.1;
gammad = 0.01;
lambdad = 0.5;

parameters.theta = theta;
parameters.gammac = gammac;
parameters.gammad = gammad;
parameters.lambdac = lambdac;
parameters.lambdad = lambdad;
parameters.Js = Js;
parameters.Jw = Jw;
parameters.delta = delta;
parameters.Mt = Mt;
parameters.Kp = Kp;
parameters.Kd = Kd;
parameters.Ki = Ki;
parameters.Omega_max = Omega_max;
parameters.taus = taus;
parameters.u = u;


%% simulate
sys = hybridPE(parameters);
sysPID = PID(parameters);
tspan = [0 70000];
jspan = [0 1000];
config = HybridSolverConfig('RelTol', 1e-10, 'AbsTol', 1e-10,  ...
                            'Refine', 10);

%% simulate
x0 = 0;
xdot0 = 0;
omega0 = 0;
tau0 = 0;
thetahat0 = 0;
psi0 = zeros(4,1);
eta0 = -[x0;xdot0;omega0;tau0];
z0 = [x0;xdot0;omega0;tau0;thetahat0;psi0;eta0];

sol = sys.solve(z0, tspan, jspan, config);

xerrInt0 = 0;
z0PID = [x0;xdot0;omega0;tau0;xerrInt0];
solPID = sysPID.solve(z0PID, tspan, jspan, config);


%%
HybridPlotBuilder.defaults.reset()

% journal
HybridPlotBuilder.defaults.set('flow line width', 3, ...
                               'jump line width', 2,...
                               'label size', 33,...
                               'tick label size', 30,...
                               't_label', '$t$ [s]')
legendFontSize = 26;

% dissertation
% HybridPlotBuilder.defaults.set('flow line width', 3, ...
%                                'jump line width', 2,...
%                                'label size', 23,...
%                                'tick label size', 23,...
%                                't_label', '$t$ [s]')   
% legendFontSize = 22;

%%
figure; clf;
hpb = HybridPlotBuilder();
hpb.flowColor('blue')...
   .labels('$| \theta - \hat \theta |$ [N-m]')...
   .plotFlows(sol.select(5),@(x) norm(theta - x))
grid on; box on;
xticks([0:8]*10^4)
yticks([0:2:6]*10^(-3))

pos = get(gcf, 'Position');
set(gcf, 'Position',  [pos(1), pos(2), 1.8*pos(3), pos(4)])
set(gca, 'LooseInset', get(gca,'TightInset'))
movegui(gcf,'center')

%% 
figure; clf;
hpb = HybridPlotBuilder();
tlt = tiledlayout(3,1,'TileSpacing','compact');
bgAx = axes(tlt,'XTick',[],'YTick',[],'Box','on');
bgAx.Layout.TileSpan = [2 1];

ax1 = axes(tlt);
ax1.Layout.Tile = 1;
hold on
% hpb.plotFlows(sol.select(1),@(x) rad2deg(u - x))
hpb.flowColor('green')...
   .tLabel('')...
   .plotFlows(solPID.select(1),@(x) rad2deg(u - x))
hpb.flowColor('blue')...
   .tLabel('')...
   .plotFlows(sol.select(1),@(x) rad2deg(u - x))
yline(ax1,1,':');
ax1.XAxis.Visible = 'off';
ax1.Box = 'off';
ax1.XAxis.Exponent = 4;
ylim(ax1,[1 5])
yticks([1 3 5])
xlim(tspan)
grid on;

ax2 = axes(tlt);
ax2.Layout.Tile = 2;
hold on
hpb.plotFlows(sol.select(1),@(x) rad2deg(u - x))
hpb.flowColor('green')...
   .legend({'PID control'},'Fontsize',legendFontSize,'Location','SouthEast','Orientation','Horizontal')...
   .plotFlows(solPID.select(1),@(x) rad2deg(u - x))
hpb.flowColor('blue')...
   .legend({'Hybrid control'},'Fontsize',legendFontSize,'Location','SouthEast','Orientation','Horizontal')...
   .plotFlows(sol.select(1),@(x) rad2deg(u - x))
hpb.reorderLegendEntries([2 1])
yline(ax2,0.0002,':');
ax2.Box = 'off';
ax2.YAxis.Exponent = 0;
% ax2.XAxis.Exponent = 4;
% ylabel('$x_{\rm des} - x$ [deg]','interpreter','latex','FontWeight','bold','FontSize',32);
% xlabel('$t$ [s]','interpreter','latex','FontWeight','bold','FontSize',26)
ylim(ax2,[-0.04 0.005])
% yticks([-0.015 -0.005 0 0.005])
xlim(tspan)
xticklabels({})
yticks([-0.04 -0.02 0])
yticklabels({'-0.04','-0.02','0'})
grid on;

ax4 = axes(tlt);
ax4.Layout.Tile = 3;
hold on
% plotflowsColorMod(t,j,convangvel(z(:,3),'rad/s','rpm'),[],[],{'color',Q(1,:)},[],3,1.5);
hpb = HybridPlotBuilder();
hpb.flowColor('green')...
   .plotFlows(solPID.select(3),@(x) convangvel(x,'rad/s','rpm'))
hpb.flowColor('blue')...
   .plotFlows(sol.select(3),@(x) convangvel(x,'rad/s','rpm'))
% ax4.XAxis.Exponent = 3;
% ylim(ax4,[96.75 100.25])
ax4.XAxis.Exponent = 4;
xlim(tspan)
ylabel('$\Omega$ [RPM]','interpreter','latex','FontWeight','bold','FontSize',32)
xlabel('$t$ [s]','interpreter','latex','FontWeight','bold','FontSize',32)
grid on; box on

linkaxes([ax1 ax2 ax4], 'x')

tlt.Padding = "none";
tlt.TileSpacing = "compact";
pos = get(gcf, 'Position');
set(gcf, 'Position',  [pos(1), pos(2), pos(3)*1.8, pos(4)*1.8])
set(bgAx.YLabel,'String','$z_{\rm des} - z$ [deg]','interpreter','latex','FontWeight','bold','FontSize',32)
bgAx.YLabel.Position(1) = bgAx.YLabel.Position(1) - 0.13;
movegui(gcf,'center')
