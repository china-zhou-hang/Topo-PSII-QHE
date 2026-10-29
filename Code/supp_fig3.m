clear all;clc;
load('fig2data/cmap_RdYlBu.mat');
cmap = cmap_RdYlBu;
set(groot,'DefaultFigureColormap',cmap);

Vlist = linspace(0.8,1.0,1000);
[P43,PEP] = fun01(Vlist);
save('supp_fig3_data.mat','Vlist','PEP');
%%
close all
color1=[218,243,241]/255;
color2=[253,235,208]/255;
figure('Position', [100,200,700,400]);
ax1= subplot(1,1,1);
patch(ax1, [0.9 1.1 1.1 0.9], [-10 -10 300 300],color1, 'FaceAlpha',1, 'EdgeColor','none');hold on;
patch(ax1, [1.05 1.0 1.0 1.05], [-10 -10 300 300],color2, 'FaceAlpha',1, 'EdgeColor','none');hold on;
plot(Vlist,PEP,'LineWidth',6, 'Color',[238,90,57]/255);hold on;
plot(0.95*ones(1,20), linspace(-0.5,3.95,20),'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(linspace(0.6,0.95,20), 3.95*ones(1,20),'--','LineWidth',2, 'Color',[150,150,150]/255);hold on;
plot(0.95, 3.95, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [238,90,57]/255, 'MarkerFaceColor', [255,255,255]/255,'LineWidth', 1.8);
ax1.Box      = 'on';
ax1.Layer    = 'top';
ax1.LineWidth = 2;
ax1.TickDir  = 'in';
set(ax1, 'FontSize', 20);
xlim([0.8 1.05])
xticks([0.8 0.9 0.95 1.0 1.05]);
xticklabels({'0.8','0.9','0.95','1.0','1.05'});
ylim([-6 300])
yticks([0 150 300]);
% yticklabels({'0','3.0','6.0'});
xlabel('$V$ (units of $\omega_{21}$)', 'Interpreter', 'latex', 'FontSize', 24); 
ylabel('$\Omega_{\rm EP}$/$\Omega_{\rm c}$','Interpreter','latex', 'FontSize', 24); 
text(0.45,0.15,'(0.95,3.95)', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])


function [P43,PEP] = fun01(Vlist)
rh = 2.23e-7*8065.543887;
rc = 2.23e-3*8065.543887;
Gc = 0.0248*8065.543887;
nh = 10000;
nc=0;
Nc=0;
Eac=1611;
Eab=14856;
n=length(Vlist);
eta=zeros(1,n);
jj=zeros(1,n);
P43=zeros(1,n);
PEP=zeros(1,n);
Geff=(rh*(nh+1) + rc*(nc+1));
for i =1:n
    V=Vlist(i);
    G = Gc*exp((1.44*8065.543887)*(1-V)/208.5);
    OmeT = sqrt(0.5 * (G + rc*nc)*0.5 * (rh*(nh+1) + rc*(nc+1)));
    c3=(G+rc+nh*rh)*(G*Gc*(rc+2*nh*rh)+(G+Gc)*nh*rc*rh);
    c4=G*(rc+rh*nh)*(G*Gc+G*nh*rh+3*Gc*nh*rh);
    Omec=sqrt(c3/c4)*OmeT;
    OmeEP=abs(Geff-G)/4;
    % PEP(1,i)=Omec/OmeEP;
    PEP(1,i)=OmeEP/Omec;
end
end
