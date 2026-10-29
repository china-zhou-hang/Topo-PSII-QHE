clear; close all; clc;
load('fig4data/fig4_data.mat');
load('fig2data/cmap_RdYlBu.mat');
cmap = cmap_RdYlBu;
set(groot,'DefaultFigureColormap',cmap);

%% ---------------- plot ----------------
color1=[141 200 255]/255;%NR1
color2=[255 169 116]/255;%NR2
color3=[0 111 188]/255;%num
color4=[255 64 36]/255;%ana1
color5=color4;%ana2
color6=[0 111 188]/255;%line3 a-b
color7=[0 0 0]/255;%;line1 互易
color8=[238 90 57]/255;%line2 b-a
color9=[150 150 150]/255;
color10=color1/1.2;
faceA=0.2;
h = figure;
set(h, 'Position', [200, 0, 1000, 1000*0.8]);

lft   = 0.11;
hgap  = 0.15;
width = 0.34;
top_bot = 0.6;
top_h   = 0.31;
bot_bot = 0.15;
xL = lft;
xR = lft + width + hgap;

ax1 = axes('Position', [xL, top_bot, width, top_h]);hold(ax1,'on');
% box(ax1,'on');
imagesc(xOmega,rcFactorList,etaMap);
set(gca, 'YDir', 'normal');
axis tight;
% contour(ax1,xOmega,rcFactorList,nuMap,[0.5 0.5],'k--','LineWidth',2.0);
% plot(sqrt(p0.G*(p0.rh*(p0.nh+1)+rcFactorList))/2/p0.omega43,rcFactorList,'--','LineWidth',2.0, 'Color',[150,150,150]/255);hold on;
plot(sqrt(p0.G*(p0.rh*(p0.nh+1)+p0.rc*rcFactorList))/2/p0.omega43,rcFactorList,'--','LineWidth',2.0, 'Color',[150,150,150]/255);hold on;
% plot(ax1,OmegaT0/p0.omega43,1,'kp','MarkerFaceColor','w','MarkerSize',10,'LineWidth',1.5);
cb1 = colorbar(ax1,'Location','northoutside');
axpos = ax1.Position;               % 主轴位置
cb1.Position(1) = axpos(1);          % 左对齐
cb1.Position(3) = axpos(3)/2;          % 宽度一致
cb1.Position(2) = axpos(2) + axpos(4) + 0.01;  % 向上偏一点
set(ax1, 'FontSize', 18);
clim(ax1, [0.02 0.7]); 
cb1.Ticks = [0.02 0.36 0.7];
cb1.TickLabels = {'', '', ''};
xlim([0 0.6])
xticks([0 0.3 0.6])
xticklabels({'0','0.3','0.6'});
ylim([0.1 1.5])
ax1.TickLabelInterpreter = 'latex';
set(ax1, 'FontSize', 18);
xlabel(ax1,'$\Omega/\omega_{32}$','Interpreter','latex', 'FontSize', 22);
ylabel(ax1,'$r_c/r_c^{(0)}$','Interpreter','latex', 'FontSize', 22);
ax1.Box      = 'on';
ax1.Layer    = 'top';
ax1.LineWidth = 2;
ax1.TickDir  = 'in';
text(0.85,0.90,'(a)', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])
% text(0.55, 1.07, 'Efficiency $\eta$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 19);
text(-0.00, 1.16, '$2\%$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.2, 1.16, '$36\%$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.43, 1.16, '$70\%$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.52, 1.06, '$\eta$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);

ax2 = axes('Position', [xR, top_bot, width, top_h]);hold(ax2,'on');
imagesc(xOmega,rcFactorList,sensitivityMap);
set(gca, 'YDir', 'normal');
axis tight;
% plot(sqrt(p0.G*(p0.rh*(p0.nh+1)+rcFactorList))/2/p0.omega43,rcFactorList,'--','LineWidth',2.0, 'Color',[150,150,150]/255);hold on;
plot(sqrt(p0.G*(p0.rh*(p0.nh+1)+p0.rc*rcFactorList))/2/p0.omega43,rcFactorList,'--','LineWidth',2.0, 'Color',[150,150,150]/255);hold on;
% contour(ax2,xOmega,rcFactorList,nuMap,[0.5 0.5],'k--','LineWidth',2.0);
% plot(ax2,OmegaT0/p0.omega43,1,'kp','MarkerFaceColor','w','MarkerSize',10,'LineWidth',1.5);
cb2 = colorbar(ax2,'Location','northoutside');
axpos = ax2.Position;               % 主轴位置
cb2.Position(1) = axpos(1);          % 左对齐
cb2.Position(3) = axpos(3)/2;          % 宽度一致
cb2.Position(2) = axpos(2) + axpos(4) + 0.01;  % 向上偏一点
set(ax2, 'FontSize', 18);
clim(ax2, [0 0.8]); 
cb2.Ticks = [0 0.4 0.8];
cb2.TickLabels = {'', '', ''};
xlim([0 0.6])
xticks([0 0.3 0.6])
xticklabels({'0','0.3','0.6'});
ylim([0.1 1.5])
ax2.TickLabelInterpreter = 'latex';
set(ax2, 'FontSize', 18);
xlabel(ax2,'$\Omega/\omega_{32}$','Interpreter','latex', 'FontSize', 22);
ylabel(ax2,'$r_c/r_c^{(0)}$','Interpreter','latex', 'FontSize', 22);
ax2.Box      = 'on';
ax2.Layer    = 'top';
ax2.LineWidth = 2;
ax2.TickDir  = 'in';
text(0.85,0.90,'(b)', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])
% text(0.55, 1.07, 'Efficiency $\eta$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 19);
text(-0.00, 1.16, '$0$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.2, 1.16, '$0.4$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
text(0.43, 1.16, '$0.8$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 16);
% text(0.52, 1.06, '$\left|\partial\ln\eta/\partial\ln r_c\right|$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
text(0.53, 1.11, '$\left|\frac{\partial\ln(\eta)}{\partial\ln(r_c)}\right|$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);

ax3 = axes('Position', [xL, bot_bot, width, top_h]);hold(ax3,'on');
h1 = plot(ax3,epsilonList,etaWorst,'LineWidth',4, 'Color',[100,136,192]/255);
h2 = plot(ax3,epsilonList,JWorst,'--','LineWidth',3.5, 'Color',[100,136,192]/255);
ylabel(ax3,'Worst-case retention','Interpreter','latex');
ylimLeft = [max(0,min([etaWorst,JWorst])-0.04),1.025];
ylim(ax3,ylimLeft);

yyaxis(ax3,'right');
h3 = plot(ax3,epsilonList,marginWorst,'-','LineWidth',4, 'Color',[238,90,57]/255);
yline(ax3,0,'--','LineWidth',2.0,'HandleVisibility','off', 'Color',[150,150,150]/255);
h4 = xline(ax3,epsilonT,'--','LineWidth',2.0, 'Color',[150,150,150]/255);
ax3.TickLabelInterpreter = 'latex';
xticks([0 0.5 1 1.5])
xticklabels({'0','','1.0','1.5'})
ylim([-0.14 0.1])
yticks([-0.1 0 0.1])
% yticklabels({'$-0.1$','$0$','$0.1$'})
set(ax3, 'FontSize', 18);
xlabel('Uncertainty radius $\epsilon$','Interpreter', 'latex', 'FontSize', 22);
% ylabel('$m^{wc}_{\rm top}/\Omega_T^{(0)}$','Interpreter', 'latex', 'FontSize', 22);
text(1.14, 0.3,'$m^{\rm wc}_{\rm top}/\Omega_T^{(0)}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24,'Color',[0,0,0]/255, 'Rotation', 90)

ax3.Box      = 'on';
ax3.Layer    = 'top';
ax3.LineWidth = 2;
ax3.TickDir  = 'in';
ax3.YColor = [0,0,0]/255;
text(0.85,0.90,'(c)', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])
text(0.27, -0.09, '$\epsilon_T$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
lgd=legend([h1,h2,h3], {'$R_{\eta}^{\rm wc}$','$R_{J}^{\rm wc}$','$m^{\rm wc}_{\rm top}/\Omega_T^{(0)}$'}, 'Interpreter', 'latex');
lgd.Position = [0.15, 0.17, 0.1, 0.1];  
lgd.Box = 'off';
lgd.FontSize = 18;
lgd.NumColumns = 1;

ax4 = axes('Position', [xR, bot_bot, width, top_h]);hold(ax4,'on');
yyaxis(ax4,'left');
h5 = fill(ax4,[xOmegaMC,fliplr(xOmegaMC)],[ratioP05,fliplr(ratioP95)],[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.60);
h6 = plot(ax4,xOmegaMC,ratioMedian,'LineWidth',3, 'Color',[100,136,192]/255);
% h7 = plot(ax4,xOmegaMC,ratioP05,'--','LineWidth',2, 'Color',[100,136,192]/255);
yline(ax4,1,':','LineWidth',1.0,'HandleVisibility','off');
% ylabel('$R_{\eta}=\eta_{\rm dis}/\eta_0$','Interpreter', 'latex', 'FontSize', 22);
text(-0.12, 0.2,'$R_{\eta}=\eta_{\rm dis}/\eta_0$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24,'Color',[0,0,0]/255, 'Rotation', 90)
ylim(ax4,[0.85 1.15]);
yticks([0.85 1.0 1.15])
ax4.YColor = [0,0,0]/255;


yyaxis(ax4,'right');
h8 = plot(ax4,xOmegaMC,probNu1,'LineWidth',4, 'Color',[238,90,57]/255);

ylim(ax4,[-0.05 1.05]);
yticks([0 1])
% yticklabels({'$-0.1$','$0$','$0.1$'})
h9 = xline(ax4,OmegaT0/p0.omega43,'--','LineWidth',2.0, 'Color',[150,150,150]/255);
box(ax4,'on');
xlabel('$\Omega/\omega_{32}$','Interpreter', 'latex', 'FontSize', 22);
% ylabel('$\Pr(\nu_{\rm eff}=1)$','Interpreter', 'latex', 'FontSize', 22);
text(1.08, 0.2,'$\Pr(\nu_{\rm eff}=1)$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24,'Color',[0,0,0]/255, 'Rotation', 90)

text(0.85,0.90,'(d)', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])
ax4.Box      = 'on';
ax4.Layer    = 'top';
ax4.LineWidth = 2;
ax4.TickDir  = 'in';
ax4.TickLabelInterpreter = 'latex';
set(ax4, 'FontSize', 18);
ax4.YColor = [0,0,0]/255;
text(0.14, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);

lgd2=legend([h5,h6,h8], {'5--95\% interval','median $R_{\eta}$','$\Pr(\nu_{\rm eff}=1)$'}, 'Interpreter', 'latex');
lgd2.Position = [0.785, 0.158, 0.1, 0.1];  
lgd2.Box = 'off';
lgd2.FontSize = 18;
lgd2.NumColumns = 1;




