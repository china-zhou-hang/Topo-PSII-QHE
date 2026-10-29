clear;clc;
load('fig2data/cmap_RdYlBu.mat')%color:-RdYlBu
cmap = cmap_RdYlBu;
set(groot,'DefaultFigureColormap',cmap);  % 之后新建的所有 figure 默认用此色表
eV2cm=8065.543887;
E0 = 0; E1 = 0.18*eV2cm; E2 = 1.62*eV2cm; E3 = 1.82*eV2cm;
H0 = diag([E0,E1,E2,E3]);
E30=E3-E0;
E21=1.44*eV2cm;
E32=0.20*eV2cm;
nh=10000; 
rh = 2.23e-7*eV2cm; 
rc=2.23e-3*eV2cm; 
Gc=0.0248*eV2cm;
kB=0.6950;
Tc=300;
Th=6000;
etaC=1-Tc/Th;
kBTc=kB*Tc;
nc=0; 
Nc=0; 
Glist = logspace(log10(1e1), log10(1e6), 100);
Geff = (rh*(nh+1) + rc);
Ocrit = sqrt(Glist .* Geff)/2;
Xcrit = Ocrit / E32;
Ycrit = (kBTc*log(Gc./Glist) + E21) / E21;
[~,ord] = sort(Xcrit);
Xcrit = Xcrit(ord);
Ycrit = Ycrit(ord);

Glist   = logspace(log10(1e1), log10(1e6), 100);
Omelist = linspace(0, 0.6, 100)*E32;
Vlist = (kBTc * log(Gc./Glist) + E21)/E21;
[X, Y] = meshgrid(Omelist/E32, Vlist);
Omelist1 = linspace(0, 0.6, 10)*E32;

V=0.95*E21;
G = Gc*exp(E21*(1-V/E21)/kBTc);
OmeT = sqrt(G*((rh*(nh+1) + rc)))/2/E32;

load('fig2data/fig2_data.mat')




etamin = G*Gc*V*rc/(E30*etaC*(G*Gc*(rc+2*nh*rh)+(G+Gc)*nh*rc*rh));
etamax = G*Gc*V/(E30*etaC*(G*(Gc+nh*rh)+3*nh*rh*Gc));
%%
% clc
close all
brcolor=[237,249,249]/255;
excolor=[254,245,232]/255;
FaceA=1;
textbrcolor=[0 0 0];
textexcolor=[0 0 0];
lft=0.1;
hgap = 0.1;  % 水平间隙
width = 0.213;  % 子图宽度
height = 0.33; % 子图高度
bot_bot=0.15;
top_bot=0.54;

h = figure('Position', [200,0,1600,900]);
ax1 = axes('Position', [lft, top_bot, width, height]); 
contourf(X, Y, Exact_eta_2D, 100, 'LineColor', 'none');hold on;
plot(Xcrit, Ycrit, 'k--', 'LineWidth', 2.5,'color',[0,0,0]/255);hold on;
xlabel('');
ylabel('$V$ (units of $\omega_{21}$)', 'Interpreter', 'latex', 'FontSize', 24); 
cb1 = colorbar(ax1,'Location','northoutside');
axpos = ax1.Position;               % 主轴位置
cb1.Position(1) = axpos(1);          % 左对齐
cb1.Position(3) = axpos(3)/2;          % 宽度一致
cb1.Position(2) = axpos(2) + axpos(4) + 0.01;  % 向上偏一点
set(ax1, 'FontSize', 18);
clim(ax1, [0.12 0.72]); 
cb1.Ticks = [0.12 0.42 0.72];
cb1.TickLabels = {'', '', ''};
ax1.LineWidth = 2;
xlim([0 0.6])
xticks([0 OmeT 0.3 0.6]);
xticklabels({''});
ylim([0.85 1.05])
yticks([0.85 0.95 1.05]);
yticklabels({'0.85','0.95','1.05'});
text(0.05,0.94,'i', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])
text(0.52, 1.06, '$\eta$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
text(-0.02, 1.15, '$12\%$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.22, 1.15, '$42\%$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.45, 1.15, '$72\%$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
plot(ax1, OmeT*ones(1,20), linspace(0.8,0.95,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(ax1, linspace(0,OmeT,20), 0.95*ones(1,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(OmeT, 0.95, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [100,136,192]/255, 'MarkerFaceColor', [255,255,255]/255,'LineWidth', 1.8);


ax2 = axes('Position', [lft+width+hgap, top_bot, width, height]); 
fill([100*ones(size(Vlist)) fliplr(Xcrit)],[Vlist flip(Vlist)],  excolor, 'FaceAlpha', 1, 'EdgeColor', 'none');hold on;
fill([0*ones(size(Vlist)) Xcrit], [flip(Vlist) Vlist], brcolor, 'FaceAlpha', 1, 'EdgeColor', 'none'); hold on;
plot(Xcrit, Ycrit, 'k--', 'LineWidth', 2.5,'color',[0,0,0]/255);hold on;
xlabel('');
ylabel('$V$ (units of $\omega_{21}$)', 'Interpreter', 'latex', 'FontSize', 24); 
set(ax2, 'FontSize', 18);
clim(ax2, [-20 110]); 
ax2.LineWidth = 2;
xlim([0 0.6])
xticks([0 OmeT 0.3 0.6]);
xticklabels({''});
ylim([0.85 1.05])
yticks([0.85 0.95 1.05]);
yticklabels({'0.85','0.95','1.05'});
plot(ax2, OmeT*ones(1,20), linspace(0.8,0.95,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(ax2, linspace(0.01,OmeT,20), 0.95*ones(1,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(OmeT, 0.95, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [100,136,192]/255, 'MarkerFaceColor', [255,255,255]/255,'LineWidth', 1.8);
text(0.05,0.94,'i', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])
text(0.1, 0.1, '$\nu=0$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28);
text(0.05, 0.25, 'Trivial', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textbrcolor);
text(0.5, 0.5, '$\nu=1$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28);
text(0.4, 0.7, 'Nontrivial', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textexcolor);


ax3 = axes('Position', [lft+2*width+2*hgap, top_bot, width, height]); 
contourf(X, Y, Cl1_2D, 100, 'LineColor', 'none');hold on;
plot(Xcrit, Ycrit, 'k--', 'LineWidth', 2.5,'color',[0,0,0]/255);hold on;
xlabel('');
ylabel('$V$ (units of $\omega_{21}$)', 'Interpreter', 'latex', 'FontSize', 24); hold on;
cb3 = colorbar(ax3,'Location','northoutside');
axpos = ax3.Position;               % 主轴位置
cb3.Position(1) = axpos(1);          % 左对齐
cb3.Position(3) = axpos(3)/2;          % 宽度一致
cb3.Position(2) = axpos(2) + axpos(4) + 0.01;  % 向上偏一点
set(ax3, 'FontSize', 18);
clim(ax3, [0 0.12]); 
cb3.Ticks = [0 0.06 0.12];
cb3.TickLabels = {'', '', ''};
ax3.LineWidth = 2;
xlim([0 0.6])
xticks([0 OmeT 0.3 0.6]);
xticklabels({''});
ylim([0.85 1.05])
yticks([0.85 0.95 1.05]);
yticklabels({'0.85','0.95','1.05'});
text(0.05,0.94,'i', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])
text(0.52, 1.06, 'Coherence $\mathcal{C}_{\ell_1}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
text(0.0, 1.15, '$0$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.2, 1.15, '$0.06$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.43, 1.15, '$0.12$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
plot(ax3, OmeT*ones(1,20), linspace(0.8,0.95,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(ax3, linspace(0,OmeT,20), 0.95*ones(1,20),'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(OmeT, 0.95, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [100,136,192]/255, 'MarkerFaceColor', [255,255,255]/255,'LineWidth', 1.8);

ax4 = axes('Position', [lft, bot_bot, width, height]); 
hold(ax4,'on');
yyaxis left
patch(ax4, [0 OmeT OmeT 0], [0 0 1 1]*100,brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');
patch(ax4, [OmeT 0.6 0.6 OmeT], [0 0 1 1]*100,excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');
plot(ax4, OmeT*ones(1,20), linspace(0,100,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);
plot(ax4, linspace(0,0.6,20),100*etamin*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);
plot(ax4, linspace(0,0.6,20),100*etamax*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);
h1=plot(Omelist/E32, 100*RWA_eta_1D, '-', 'DisplayName','Non-RWA','LineWidth',6, 'Color',[100,136,192]/255);
h2=plot(Omelist1/E32, 100*Exact_eta_1D,  'o','LineWidth',1.8,'MarkerSize', 15, 'Color',[16,166,122]/255,'MarkerFaceColor', [100,136,192]/255);
ax4 = gca;
ax4.YColor = [0 0 0];
ylim(ax4, [20 80]);
yticks(ax4, [20 100*etamin 50 100*etamax 80]);
yticklabels(ax4, {'20%','','50%','','80%'});
text(-0.23, 0.48,'$\eta$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)

yyaxis right
h3=plot(Omelist/E32, 1e5*RWA_detadO2_1D, '-', 'DisplayName','Non-RWA','LineWidth',6, 'Color',[238,90,57]/255);
text(1.08, 0.4,'$\mathcal{H}_{\eta}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24,'Color',[0,0,0]/255, 'Rotation', 90)
ylim(ax4, [-1 3]);
yticks(ax4, [-1 0 3]);
yticklabels(ax4, {'-1','0','3'});
ax = gca;
ax.YColor = [0,0,0]/255;
text(1.085, 0.55,'($10^{-5}$)','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)

ax4.Box      = 'on';
ax4.Layer    = 'top';
ax4.LineWidth = 2;
ax4.TickDir  = 'in';
set(ax4, 'FontSize', 18);
xlim(ax4, [0 0.6]);
xticks([0 OmeT 0.3 0.6]);
xticklabels({'0','','0.3','0.6'});
text(ax4, 0.05, 0.94, 'ii', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);

text(-0.2, 0.14, '$\eta_{\mathrm{min}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
text(-0.2, 0.88, '$\eta_{\mathrm{max}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
text(0.28, 0.92, 'V=0.95$\omega_{21}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 20);
text(0.14, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
text(0.23, 0.4, 'TPT', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 18, 'Rotation', -90);
lgd=legend([h1, h2,h3], {'RWA','Exact','$\mathcal{H}_{\eta}$'}, 'Interpreter', 'latex');
lgd.Position = [0.20, 0.27, 0.1, 0.1];  
lgd.Box = 'off';
lgd.FontSize = 22;


ax5 = axes('Position', [lft+width+hgap, bot_bot, width, height]); 
hold(ax5,'on');
nu = ones(size(Omelist));
nu(Omelist/E32 < OmeT) = 0;
nu(Omelist/E32 >=OmeT) = 1; 
patch(ax5, [0 OmeT OmeT 0], [-1 -1 2 2],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');
patch(ax5, [OmeT 0.6 0.6 OmeT], [-1 -1 2 2],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');
plot(ax5, OmeT*ones(1,20), linspace(-1,2,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);
h2=plot(Omelist/E32,nu, 'Color',[100,136,192]/255,'LineWidth', 5);
ax5.Box      = 'on';
ax5.Layer    = 'top';
ax5.LineWidth = 2;
ax5.TickDir  = 'in';
set(ax5, 'FontSize', 18);
xlim(ax5, [0 0.6]);
ylim(ax5, [-1 2]);
xticks([0 OmeT 0.3 0.6]);
xticklabels({'0','','0.3','0.6'});
text(ax5, 0.05, 0.94, 'ii', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);
text(-0.2, 0,'Winding number $\nu$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)
text(0.28, 0.92, 'V=0.95$\omega_{21}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 20);
text(0.14, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
text(0.23, 0.4, 'TPT', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 18, 'Rotation', -90);
set([h1 h2],'HandleVisibility','off');

ax6 = axes('Position', [lft+2*width+2*hgap, bot_bot, width, height]); 
hold(ax6,'on');
yyaxis left
patch(ax6, [0 OmeT OmeT 0], [-0.01 -0.01 0.05 0.05],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');
patch(ax6, [OmeT 0.6 0.6 OmeT], [-0.01 -0.01 0.05 0.05],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');
plot(ax6, OmeT*ones(1,20), linspace(-0.01,0.05,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);
h1=plot(Omelist/E32, RWA_Cl1_1D, '-', 'DisplayName','Non-RWA','LineWidth',6, 'Color',[100,136,192]/255);
ylim(ax6, [-0.002 0.042]);
yticks(ax6, [0 0.02 0.04]);
yticklabels(ax6, {'0', '0.02','0.04'});
ax = gca;
ax.YColor = [0 0 0];
text(-0.2, 0.15,'Coherence $\mathcal{C}_{\ell_1}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)

yyaxis right
h2=plot(Omelist/E32, RWA_Jc_1D, '-', 'DisplayName','Non-RWA','LineWidth',6, 'Color',[238,90,57]/255);
ylim(ax6, [-0.8 16]);
yticks(ax6, [0 8 16]);
yticklabels(ax6, {'0', '','16'});
ax.YColor = [0,0,0]/255;
text(1.08, 0.45,'$J_{c}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)

ax6.Box      = 'on';
ax6.Layer    = 'top';
ax6.LineWidth = 2;
ax6.TickDir  = 'in';
set(ax6, 'FontSize', 18);
xlim(ax6, [0 0.6]);
xticks([0 OmeT 0.3 0.6]);
xticklabels({'0','','0.3','0.6'});
text(ax6, 0.05, 0.94, 'ii', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);
text(0.28, 0.92, 'V=0.95$\omega_{21}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 20);
text(0.14, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
text(0.23, 0.4, 'TPT', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 18, 'Rotation', -90);
lgd=legend([h1, h2], {'$\mathcal{C}_{\ell_1}$','$J_{\rm c}$'}, 'Interpreter', 'latex');
lgd.Position = [0.8, 0.2, 0.1, 0.1];  
lgd.Box = 'off';
lgd.FontSize = 22;