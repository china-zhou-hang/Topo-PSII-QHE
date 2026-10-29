clear all;clc;
load('fig2data/cmap_RdYlBu.mat');
cmap = cmap_RdYlBu;
set(groot,'DefaultFigureColormap',cmap);
nh=10000;
rh = 2.23e-7*8065.543887;
rc = 2.23e-3*8065.543887;
Gc= 0.0248*8065.543887;
E41=1.82*8065.543887;
etaC=0.95;
tic
V1=0.94;
V2=0.96;
V3=0.98;
V4=1.00;
G1 = Gc*exp((1.44*8065.543887)*(1-V1)/208.5);
G2 = Gc*exp((1.44*8065.543887)*(1-V2)/208.5);
G3 = Gc*exp((1.44*8065.543887)*(1-V3)/208.5);
G4 = Gc*exp((1.44*8065.543887)*(1-V4)/208.5);
OmeT1 = sqrt(0.5*G1*(0.5*(rh*(nh+1) + rc)))/(0.20*8065.543887);
OmeT2 = sqrt(0.5*G2*(0.5*(rh*(nh+1) + rc)))/(0.20*8065.543887);
OmeT3 = sqrt(0.5*G3*(0.5*(rh*(nh+1) + rc)))/(0.20*8065.543887);
OmeT4 = sqrt(0.5*G4*(0.5*(rh*(nh+1) + rc)))/(0.20*8065.543887);
[etalist1,rho43list1,Jcohlist1,detadO2list1] = fun01(V1,nh,rh,rc,Gc,E41,etaC);
[etalist2,rho43list2,Jcohlist2,detadO2list2] = fun01(V2,nh,rh,rc,Gc,E41,etaC);
[etalist3,rho43list3,Jcohlist3,detadO2list3] = fun01(V3,nh,rh,rc,Gc,E41,etaC);
[etalist4,rho43list4,Jcohlist4,detadO2list4] = fun01(V4,nh,rh,rc,Gc,E41,etaC);
% Omelist  = logspace(-2, log10(30),100);
Omelist  = linspace(0, 0.6,100);
etamin1 = G1*Gc*(1.44*8065.543887)*V1*rc/(E41*etaC*(G1*Gc*(rc+2*nh*rh)+(G1+Gc)*nh*rc*rh));
etamax1 = G1*Gc*(1.44*8065.543887)*V1/(E41*etaC*(G1*(Gc+nh*rh)+3*nh*rh*Gc));

etamin2 = G2*Gc*(1.44*8065.543887)*V2*rc/(E41*etaC*(G2*Gc*(rc+2*nh*rh)+(G2+Gc)*nh*rc*rh));
etamax2 = G2*Gc*(1.44*8065.543887)*V2/(E41*etaC*(G2*(Gc+nh*rh)+3*nh*rh*Gc));

etamin3 = G3*Gc*(1.44*8065.543887)*V3*rc/(E41*etaC*(G3*Gc*(rc+2*nh*rh)+(G3+Gc)*nh*rc*rh));
etamax3 = G3*Gc*(1.44*8065.543887)*V3/(E41*etaC*(G3*(Gc+nh*rh)+3*nh*rh*Gc));

etamin4 = G4*Gc*(1.44*8065.543887)*V4*rc/(E41*etaC*(G4*Gc*(rc+2*nh*rh)+(G4+Gc)*nh*rc*rh));
etamax4 = G4*Gc*(1.44*8065.543887)*V4/(E41*etaC*(G4*(Gc+nh*rh)+3*nh*rh*Gc));
toc
save('supp_fig5_data.mat','Omelist','V1','V2','V3','V4','OmeT1','OmeT2','OmeT3','OmeT4','etalist1','etalist2','etalist3','etalist4','rho43list1','rho43list2','rho43list3','rho43list4','Jcohlist1','Jcohlist2','Jcohlist3','Jcohlist4','detadO2list1','detadO2list2','detadO2list3','detadO2list4','etamin1','etamin2','etamin3','etamin4','etamax1','etamax2','etamax3','etamax4');
%% 
clc
close all
brcolor=[237,249,249]/255;
excolor=[254,245,232]/255;
FaceA=1;
textbrcolor=[0 0 0];
textexcolor=[0 0 0];
lft=0.1;
hgap = 0.08;  % 水平间隙
width = 0.15;  % 子图宽度
height = 0.33; % 子图高度
bot_bot=0.15;
top_bot=0.54;

h = figure('Position', [900,150,2000,800]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1 = axes('Position', [lft, top_bot, width, height]); 
hold(ax1,'on');
yyaxis left
patch(ax1, [0 OmeT1 OmeT1 0], [20 20 80 80],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
patch(ax1, [OmeT1 0.6 0.6 OmeT1], [20 20 80 80],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
plot(ax1, OmeT1*ones(1,20), linspace(20,80,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);hold on;
plot(ax1, linspace(0,0.6,20),100*etamin1*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(ax1, linspace(0,0.6,20),100*etamax1*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
h1=plot(Omelist, 100*etalist1, '-','LineWidth',6, 'Color',[100,136,192]/255); hold on;
text(-0.25, 0.45,'$\eta$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90);
ax1 = gca;
ax1.YColor = [0 0 0];
ylim(ax1, [20 80]);
yticks(ax1, [20 100*etamin1 50 100*etamax1 80]);
yticklabels(ax1, {'20%','','50%','','80%'});
text(ax1, 0.03, 0.94, 'i', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
text(0.3, 1.08,'V=0.94$\omega_{21}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255)

yyaxis right
h2=plot(Omelist, 1e6*detadO2list1, '-','LineWidth',6, 'Color',[238,90,57]/255); 
text(1.08, 0.4,'$\mathcal{H}_{\eta}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24,'Color',[0,0,0]/255, 'Rotation', 90)
ylim(ax1, [-4 16]);
yticks(ax1, [-4 0 16]);
yticklabels(ax1, {'-4','0','16'});
text(1.085, 0.55,'($10^{-6}$)','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)
ax1 = gca;
ax1.YColor = [0 0 0];
ax1.Box      = 'on';
ax1.Layer    = 'top';
ax1.LineWidth = 2;
ax1.TickDir  = 'in';
set(ax1, 'FontSize', 18);
xlim(ax1, [0 0.6]);
xticks([0 OmeT1 0.3 0.6]);
xticklabels({''});
text(-0.2, (100*etamin1-20)/60, '$\eta_{\mathrm{min}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
text(-0.2, (100*etamax1-20)/60, '$\eta_{\mathrm{max}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
lgd1=legend([h1,h2], {'$\eta$','$\mathcal{H}_{\eta}$'}, 'Interpreter', 'latex');
lgd1.Position = [0.16, 0.65, 0.1, 0.1];  
lgd1.Box = 'off';
lgd1.FontSize = 22;
text(-0.2,1.2, '(a)', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textexcolor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ax2 = axes('Position', [lft+width+hgap, top_bot, width, height]); 
hold(ax2,'on');
yyaxis left
patch(ax2, [0 OmeT2 OmeT2 0], [20 20 80 80],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
patch(ax2, [OmeT2 0.6 0.6 OmeT2], [20 20 80 80],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
plot(ax2, OmeT2*ones(1,20), linspace(20,80,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);hold on;
plot(ax2, linspace(0,0.6,20),100*etamin2*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(ax2, linspace(0,0.6,20),100*etamax2*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
h3=plot(Omelist, 100*etalist2, '-','LineWidth',6, 'Color',[100,136,192]/255); hold on;
text(-0.25, 0.45,'$\eta$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90);
ax2 = gca;
ax2.YColor = [0 0 0];
ylim(ax2, [20 80]);
yticks(ax2, [20 100*etamin2 50 100*etamax2 80]);
yticklabels(ax2, {'20%','','50%','','80%'});
text(ax2, 0.03, 0.94, 'i', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
text(0.3, 1.08,'V=0.96$\omega_{21}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255)

yyaxis right
h4=plot(Omelist, 1e5*detadO2list2, '-','LineWidth',6, 'Color',[238,90,57]/255); 
text(1.08, 0.4,'$\mathcal{H}_{\eta}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24,'Color',[0,0,0]/255, 'Rotation', 90)
ylim(ax2, [-2 5]);
yticks(ax2, [-2 0 5]);
yticklabels(ax2, {'-2','0','5'});
text(1.085, 0.55,'($10^{-5}$)','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)
ax2 = gca;
ax2.YColor = [0 0 0];
ax2.Box      = 'on';
ax2.Layer    = 'top';
ax2.LineWidth = 2;
ax2.TickDir  = 'in';
set(ax2, 'FontSize', 18);
xlim(ax2, [0 0.6]);
xticks([0 OmeT2 0.3 0.6]);
xticklabels({''});
text(-0.2, (100*etamin2-20)/60, '$\eta_{\mathrm{min}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
text(-0.2, (100*etamax2-20)/60, '$\eta_{\mathrm{max}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
lgd2=legend([h3,h4], {'$\eta$','$\mathcal{H}_{\eta}$'}, 'Interpreter', 'latex');
lgd2.Position = [0.39, 0.65, 0.1, 0.1];  
lgd2.Box = 'off';
lgd2.FontSize = 22;
text(-0.2,1.2, '(b)', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textexcolor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ax3 = axes('Position', [lft+2*width+2*hgap, top_bot, width, height]); 
hold(ax3,'on');
yyaxis left
patch(ax3, [0 OmeT3 OmeT3 0], [20 20 80 80],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
patch(ax3, [OmeT3 0.6 0.6 OmeT3], [20 20 80 80],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
plot(ax3, OmeT3*ones(1,20), linspace(20,80,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);hold on;
plot(ax3, linspace(0,0.6,20),100*etamin3*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(ax3, linspace(0,0.6,20),100*etamax3*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
h5=plot(Omelist, 100*etalist3, '-','LineWidth',6, 'Color',[100,136,192]/255); hold on;
text(-0.25, 0.45,'$\eta$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90);
ax3 = gca;
ax3.YColor = [0 0 0];
ylim(ax3, [20 80]);
yticks(ax3, [8 20 32]);
yticks(ax3, [20 100*etamin3 50 100*etamax3 80]);
yticklabels(ax3, {'20%','','50%','','80%'});
text(ax3, 0.03, 0.94, 'i', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
text(0.3, 1.08,'V=0.98$\omega_{21}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255)


yyaxis right
h6=plot(Omelist, 1e5*detadO2list3, '-','LineWidth',6, 'Color',[238,90,57]/255); 
text(1.08, 0.4,'$\mathcal{H}_{\eta}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24,'Color',[0,0,0]/255, 'Rotation', 90)
ylim(ax3, [-5 15]);
yticks(ax3, [-5 0 15]);
yticklabels(ax3, {'-5','0','15'});
text(1.085, 0.55,'($10^{-5}$)','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)
ax3 = gca;
ax3.YColor = [0 0 0];
ax3.Box      = 'on';
ax3.Layer    = 'top';
ax3.LineWidth = 2;
ax3.TickDir  = 'in';
set(ax3, 'FontSize', 18);
xlim(ax3, [0 0.6]);
xticks([0 OmeT3 0.3 0.6]);
xticklabels({''});
text(-0.2, (100*etamin3-20)/60, '$\eta_{\mathrm{min}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
text(-0.2, (100*etamax3-20)/60, '$\eta_{\mathrm{max}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
text(-0.2,1.2, '(c)', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textexcolor);
lgd3=legend([h5,h6], {'$\eta$','$\mathcal{H}_{\eta}$'}, 'Interpreter', 'latex');
lgd3.Position = [0.62, 0.65, 0.1, 0.1];  
lgd3.Box = 'off';
lgd3.FontSize = 22;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ax4 = axes('Position', [lft+3*width+3*hgap, top_bot, width, height]); 
hold(ax4,'on');
yyaxis left
patch(ax4, [0 OmeT4 OmeT4 0], [20 20 80 80],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
patch(ax4, [OmeT4 0.6 0.6 OmeT4], [20 20 80 80],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
plot(ax4, OmeT4*ones(1,20), linspace(20,80,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);hold on;
plot(ax4, linspace(0,0.6,20),100*etamin4*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(ax4, linspace(0,0.6,20),100*etamax4*ones(1,20) ,'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
h5=plot(Omelist, 100*etalist4, '-','LineWidth',6, 'Color',[100,136,192]/255); hold on;
text(-0.25, 0.45,'$\eta$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90);
ax4 = gca;
ax4.YColor = [0 0 0];
ylim(ax4, [20 80]);
yticks(ax4, [8 20 32]);
yticks(ax4, [20 100*etamin4 50 100*etamax4 80]);
yticklabels(ax4, {'20%','','50%','','80%'});
text(ax4, 0.03, 0.94, 'i', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
text(0.3, 1.08,'V=1.0$\omega_{21}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255)


yyaxis right
h6=plot(Omelist, 1e4*detadO2list4, '-','LineWidth',6, 'Color',[238,90,57]/255); 
text(1.08, 0.4,'$\mathcal{H}_{\eta}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24,'Color',[0,0,0]/255, 'Rotation', 90)
ylim(ax4, [-1 4]);
yticks(ax4, [-1 0 4]);
yticklabels(ax4, {'-1','0','4'});
text(1.085, 0.55,'($10^{-4}$)','Interpreter','latex', 'Units', 'normalized', 'FontSize', 22,'Color',[0,0,0]/255, 'Rotation', 90)
ax4 = gca;
ax4.YColor = [0 0 0];
ax4.Box      = 'on';
ax4.Layer    = 'top';
ax4.LineWidth = 2;
ax4.TickDir  = 'in';
set(ax4, 'FontSize', 18);
xlim(ax4, [0 0.6]);
xticks([0 OmeT4 0.3 0.6]);
xticklabels({''});
text(-0.2, (100*etamin4-20)/60, '$\eta_{\mathrm{min}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
text(-0.2, (100*etamax4-20)/60, '$\eta_{\mathrm{max}}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 24);
text(-0.2,1.2, '(d)', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 28,Color=textexcolor);
lgd4=legend([h5,h6], {'$\eta$','$\mathcal{H}_{\eta}$'}, 'Interpreter', 'latex');
lgd4.Position = [0.85, 0.65, 0.1, 0.1];  
lgd4.Box = 'off';
lgd4.FontSize = 22;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ax5 = axes('Position', [lft, bot_bot, width, height]); 
hold(ax5,'on');
yyaxis left
patch(ax5, [0 OmeT1 OmeT1 0], [-0.001 -0.001 0.04 0.04],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
patch(ax5, [OmeT1 0.6 0.6 OmeT1], [-0.001 -0.001 0.04 0.04],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
plot(ax5, OmeT1*ones(1,20), linspace(-0.001,0.04,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);hold on;
h7=plot(Omelist, rho43list1, '-','LineWidth',6, 'Color',[100,136,192]/255); hold on;
ax5.YColor = [0 0 0];
ylim(ax5, [-0.002 0.04]);
yticks(ax5, [0 0.02 0.04]);
yticklabels(ax5, {'0', '2','4'});
text(-0.2, 0.45,'$\mathcal{C}_{\ell_1}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90)
text(ax5, 0.03, 0.94, 'ii', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
text(ax5, 0.03, 1.04, '$\times 10^{-2}$','Interpreter','latex', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);

yyaxis right
h8=plot(Omelist, Jcohlist1, '-','LineWidth',6, 'Color',[238,90,57]/255); hold on;
ylim(ax5, [-0.5 16]);
yticks(ax5, [0 16]);
yticklabels(ax5, {'0','16'});
ax5.YColor = [0,0,0]/255;
text(1.08, 0.45,'$J_{c}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90)
ax5 = gca;
ax5.YColor = [0 0 0];
ax5.Box      = 'on';
ax5.Layer    = 'top';
ax5.LineWidth = 2;
ax5.TickDir  = 'in';
set(ax5, 'FontSize', 18);
xlim(ax5, [0 0.6]);
xticks([0 OmeT1 0.3 0.6]);
xticklabels({'0','','0.3','0.6'});
xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);
text(OmeT1/0.6, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
lgd5=legend([h7, h8], {'$\mathcal{C}_{\ell_1}$','$J_{\rm c}$'}, 'Interpreter', 'latex');
lgd5.Position = [0.16, 0.2, 0.1, 0.1];  
lgd5.Box = 'off';
lgd5.FontSize = 22;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ax6 = axes('Position', [lft+width+hgap, bot_bot, width, height]); 
hold(ax6,'on');
yyaxis left
patch(ax6, [0 OmeT2 OmeT2 0], [-0.002 -0.002 0.14 0.14],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
patch(ax6, [OmeT2 0.6 0.6 OmeT2], [-0.002 -0.002 0.14 0.14],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
plot(ax6, OmeT2*ones(1,20), linspace(-0.002,0.14,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);hold on;
h9=plot(Omelist, rho43list2, '-','LineWidth',6, 'Color',[100,136,192]/255); hold on;
ax6.YColor = [0 0 0];
ylim(ax6, [-0.004 0.08]);
yticks(ax6, [0 0.04 0.08]);
yticklabels(ax6, {'0', '4','8'});
text(ax6, 0.03, 1.04, '$\times 10^{-2}$','Interpreter','latex', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
text(-0.2, 0.45,'$\mathcal{C}_{\ell_1}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90)
text(ax6, 0.03, 0.94, 'ii', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);

yyaxis right
h10=plot(Omelist, Jcohlist2, '-','LineWidth',6, 'Color',[238,90,57]/255); hold on;
ylim(ax6, [-0.5 16]);
yticks(ax6, [0 16]);
yticklabels(ax6, {'0','16'});
ax6.YColor = [0,0,0]/255;
text(1.08, 0.45,'$J_{c}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90)
ax6 = gca;
ax6.YColor = [0 0 0];
ax6.Box      = 'on';
ax6.Layer    = 'top';
ax6.LineWidth = 2;
ax6.TickDir  = 'in';
set(ax6, 'FontSize', 18);
xlim(ax6, [0 0.6]);
xticks([0 OmeT2 0.3 0.6]);
xticklabels({'0','','0.3','0.6'});
xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);
text(OmeT2/0.6, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
lgd5=legend([h9, h10], {'$\mathcal{C}_{\ell_1}$','$J_{\rm c}$'}, 'Interpreter', 'latex');
lgd5.Position = [0.36, 0.2, 0.1, 0.1];  
lgd5.Box = 'off';
lgd5.FontSize = 22;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ax7 = axes('Position', [lft+2*width+2*hgap, bot_bot, width, height]); 
hold(ax7,'on');
yyaxis left
patch(ax7, [0 OmeT3 OmeT3 0], [-0.003 -0.003 0.16 0.16],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
patch(ax7, [OmeT3 0.6 0.6 OmeT3], [-0.003 -0.003 0.16 0.16],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
plot(ax7, OmeT3*ones(1,20), linspace(-0.003,0.16,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);hold on;
h11=plot(Omelist, rho43list3, '-','LineWidth',6, 'Color',[100,136,192]/255); hold on;
ax7.YColor = [0 0 0];
ylim(ax7, [-0.006 0.12]);
yticks(ax7, [0 0.06 0.12]);
yticklabels(ax7, {'0', '6','12'});
text(ax7, 0.03, 1.04, '$\times 10^{-2}$','Interpreter','latex', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
text(-0.2, 0.45,'$\mathcal{C}_{\ell_1}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90)
text(ax7, 0.03, 0.94, 'ii', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);

yyaxis right
h12=plot(Omelist, Jcohlist3, '-','LineWidth',6, 'Color',[238,90,57]/255); hold on;
ylim(ax7, [-0.5 16]);
yticks(ax7, [0 16]);
yticklabels(ax7, {'0','16'});
ax7.YColor = [0,0,0]/255;
text(1.08, 0.45,'$J_{c}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90)
ax7 = gca;
ax7.YColor = [0 0 0];
ax7.Box      = 'on';
ax7.Layer    = 'top';
ax7.LineWidth = 2;
ax7.TickDir  = 'in';
set(ax7, 'FontSize', 18);
xlim(ax7, [0 0.6]);
xticks([0 OmeT3 0.3 0.6]);
xticklabels({'0','','0.3','0.6'});
xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);
text(OmeT3/0.6, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
lgd4=legend([h11, h12], {'$\mathcal{C}_{\ell_1}$','$J_{\rm c}$'}, 'Interpreter', 'latex');
lgd4.Position = [0.62, 0.28, 0.1, 0.1];  
lgd4.Box = 'off';
lgd4.FontSize = 22;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
ax8 = axes('Position', [lft+3*width+3*hgap, bot_bot, width, height]); 
hold(ax8,'on');
yyaxis left
patch(ax8, [0 OmeT4 OmeT4 0], [-0.004 -0.004 0.18 0.18],brcolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
patch(ax8, [OmeT4 0.6 0.6 OmeT4], [-0.004 -0.004 0.18 0.18],excolor, 'FaceAlpha',FaceA, 'EdgeColor','none');hold on;
plot(ax8, OmeT4*ones(1,20), linspace(-0.004,0.18,20),'--','LineWidth',2, 'Color',[100,136,192, 0.5*255]/255);hold on;
h11=plot(Omelist, rho43list4, '-','LineWidth',6, 'Color',[100,136,192]/255); hold on;
ax8.YColor = [0 0 0];
ylim(ax8, [-0.008 0.16]);
yticks(ax8, [0 0.08 0.16]);
yticklabels(ax8, {'0', '8','16'});
text(ax8, 0.03, 1.04, '$\times 10^{-2}$','Interpreter','latex', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);
text(-0.2, 0.45,'$\mathcal{C}_{\ell_1}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90)
text(ax8, 0.03, 0.94, 'ii', 'Units','normalized', 'FontSize', 24, 'Color',[0 0 0]);

yyaxis right
h12=plot(Omelist, Jcohlist4, '-','LineWidth',6, 'Color',[238,90,57]/255); hold on;
ylim(ax8, [-0.5 16]);
yticks(ax8, [0 16]);
yticklabels(ax8, {'0','16'});
ax8.YColor = [0,0,0]/255;
text(1.08, 0.45,'$J_{c}$','Interpreter','latex', 'Units', 'normalized', 'FontSize', 26,'Color',[0,0,0]/255, 'Rotation', 90)
ax8 = gca;
ax8.YColor = [0 0 0];
ax8.Box      = 'on';
ax8.Layer    = 'top';
ax8.LineWidth = 2;
ax8.TickDir  = 'in';
set(ax8, 'FontSize', 18);

xlim(ax8, [0 0.6]);
xticks([0 OmeT4 0.3 0.6]);
xticklabels({'0','','0.3','0.6'});
xlabel('$\Omega$ (units of $\omega_{32}$)','Interpreter','latex', 'FontSize', 24);
text(OmeT4/0.6, -0.09, '$\Omega_{\mathrm{T}}$', 'Interpreter', 'latex', 'Units', 'normalized','FontSize', 20);
lgd4=legend([h11, h12], {'$\mathcal{C}_{\ell_1}$','$J_{\rm c}$'}, 'Interpreter', 'latex');
lgd4.Position = [0.85, 0.28, 0.1, 0.1];  
lgd4.Box = 'off';
lgd4.FontSize = 22;



%% 
function [etalist,rho43list,Jcohlist,detadO2list] = fun01(V,nh,rh,rc,Gc,E41,etaC)
G= Gc*exp((1.44*8065.543887)*(1-V)/208.5);
Omelist  = linspace(0, 0.6,100)*(0.20*8065.543887);
nO = numel(Omelist);
rho43list = zeros(1,nO);
Jcohlist  = zeros(1,nO);
etalist   = zeros(1,nO);
detadO2list  = zeros(1,nO);
for iO = 1:nO
Omeac = Omelist(iO);
Wh=nh*rh; Rh=rh*(nh+1); kappa=(G+Rh+rc)/2; n0=G*Gc+G*Wh+Gc*Rh+2*Gc*Wh; d0=G*Gc*(Rh+Wh+rc)+Wh*rc*(G+Gc); OmegaC2=kappa*d0/(2*n0);
etamin=G*Gc*(1.44*8065.543887)*V*rc/(E41*etaC*(G*Gc*(rc+2*Wh)+(G+Gc)*Wh*rc));
etamax=G*Gc*(1.44*8065.543887)*V/(E41*etaC*(G*(Gc+Wh)+3*Wh*Gc));
rho43=1i*Gc*Wh*(G-rc)*Omeac/(kappa*d0+2*n0*Omeac^2);
etalist(iO)=etamin+(etamax-etamin)*Omeac^2/(OmegaC2+Omeac^2);
rho43list(iO)=2*abs(rho43);
Jcohlist(iO)=2*Omeac*imag(rho43);
detadO2list(iO)=2*(etamax-etamin)*OmegaC2*(OmegaC2-3*Omeac^2)/(OmegaC2+Omeac^2)^3;
end
% y=rho43list;
% rho43list = y;
% k = 30;
% y0 = y(1);
% t = linspace(0,1,k)';           % 0→1
% w = 0.5*(1+cos(pi*t));          % 1→0 的余弦窗（Hann 半窗）
% rho43list(1:k) = y(1:k) - w' * y0;      % 仅改前10点
% rho43list(1) = 0;

% y=Jcohlist;
% Jcohlist = y;
% k = 30;
% y0 = y(1);
% t = linspace(0,1,k)';           % 0→1
% w = 0.5*(1+cos(pi*t));          % 1→0 的余弦窗（Hann 半窗）
% Jcohlist(1:k) = y(1:k) - w' * y0;      % 仅改前10点
% Jcohlist(1) = 0;
end