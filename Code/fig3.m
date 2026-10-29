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
Vlist = linspace(0.8,1.0,100);
[P43,Pstar] = fun01(Vlist,nh,rh,rc,Gc,E21,kBTc);
save('fig3_data.mat','Vlist','P43','Pstar');
%%
close all
color1=[218,243,241]/255;
color2=[253,235,208]/255;
h = figure('Position', [200,200,700,600]);
ax1 = axes('Position', [0.2,0.58,0.75,0.35]);
patch(ax1, [0.9 1.1 1.1 0.9], [-0.18 -0.18 -0.14 -0.14],color1, 'FaceAlpha',1, 'EdgeColor','none');hold on;
patch(ax1, [1.05 1.0 1.0 1.05], [-0.18 -0.18 -0.14 -0.14],color2, 'FaceAlpha',1, 'EdgeColor','none');hold on;
plot(Vlist,P43,'LineWidth',6, 'Color',[100,136,192]/255);hold on;
plot(0.95*ones(1,20), linspace(-0.17,0.14,20),'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
ax1.Box      = 'on';
ax1.Layer    = 'top';
ax1.LineWidth = 2;
ax1.TickDir  = 'in';
set(ax1, 'FontSize', 20);
xlim([0.8 1.05])
xticks([0.8 0.9 0.95 1.0 1.05]);
xticklabels({''});
ylim([-0.1615 -0.14])
yticks([-0.16 -0.15 -0.14]);
yticklabels({'-0.16','-0.15','-0.14'});
ylabel('$R_{T}$','Interpreter','latex', 'FontSize', 24); 
text(0.9,0.9,'(a)', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])
plot(0.95, -0.159, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [100,136,192]/255, 'MarkerFaceColor', [255,255,255]/255,'LineWidth', 1.8);

ax2 = axes('Position', [0.2,0.15,0.75,0.35]);
patch(ax2, [0.9 1.1 1.1 0.9], [-0.02 -0.02 0.01 0.01],color1, 'FaceAlpha',1, 'EdgeColor','none');hold on;
patch(ax2, [1.05 1.0 1.0 1.05], [-0.02 -0.02 0.01 0.01],color2, 'FaceAlpha',1, 'EdgeColor','none');hold on;
plot(Vlist,Pstar,'LineWidth',6, 'Color',[238,90,57]/255);hold on;
plot(0.95*ones(1,20), linspace(-0.5,0.07,20),'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(0.95,-0.0004, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [238,90,57]/255, 'MarkerFaceColor', [255,255,255]/255,'LineWidth', 1.8);
ax2.Box      = 'on';
ax2.Layer    = 'top';
ax2.LineWidth = 2;
ax2.TickDir  = 'in';
set(ax2, 'FontSize', 20);
xlim([0.8 1.05])
xticks([0.8 0.9 0.95 1.0 1.05]);
xticklabels({'0.8','0.9','0.95','1.0','1.05'});
xlabel('$V$ (units of $\omega_{21}$)', 'Interpreter', 'latex', 'FontSize', 24); 
ylim([-0.02 0.001])
yticks([-0.02 -0.01 0]);
yticklabels({'-0.02','-0.01','0'});
ylabel('$R^{\ast}_{T}$','Interpreter','latex', 'FontSize', 24); 
text(0.9,0.9,'(b)', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])


function [P43,Pstar] = fun01(Vlist,nh,rh,rc,Gc,E21,kBTc)
n=length(Vlist);
P43=zeros(1,n);
Pstar=zeros(1,n);
Wh=nh*rh;
for i =1:n
    V=Vlist(i);
    G = Gc*exp(E21*(1-V)/kBTc);
    Geff=(rh*(nh+1) + rc);
    zeta3=(G+rc+nh*rh)*(G*Gc*(rc+2*nh*rh)+(G+Gc)*nh*rc*rh);
    zeta4=G*(rc+rh*nh)*(G*Gc+nh*rh*(G+3*Gc));
    OmeT = sqrt(G*Geff)/2;
    OmeC=sqrt(zeta3/zeta4)*OmeT;
    alpha_G=((3*Gc+Wh)*G+Wh*(Gc + Wh))/(2*(G+Wh)*(Gc + Wh));
    OmeT1 =sqrt(alpha_G*G*Geff)/2;
    P43(1,i)=(OmeT-OmeC)/(OmeC);
    Pstar(1,i)=(OmeT1-OmeC)/(OmeC);
end
end
