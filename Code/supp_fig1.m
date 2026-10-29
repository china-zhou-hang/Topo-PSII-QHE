clear all;clc;
load('fig2data/cmap_RdYlBu.mat');
cmap = cmap_RdYlBu;
set(groot,'DefaultFigureColormap',cmap);
nh = 10000;
rh = 2.23e-7*8065.543887;
rc = 2.23e-3*8065.543887;
Gc = 0.0248*8065.543887;
Vlist = linspace(0,1.2,1000);
[etac,eta0,etat,eta2t] = fun01(Vlist,nh,rh,rc,Gc);
save('supp_fig1_data.mat','Vlist','etac','eta0','etat','eta2t');
%%
close all
color1=[218,243,241]/255;
h = figure('Position', [100,200,700,400]);
ax = axes('Position', [0.15,0.2,0.65,0.7]); 
patch(ax, [0.9 1 1 0.9], [-10 -10 600 600],color1, 'FaceAlpha',1, 'EdgeColor','none');hold on;
yyaxis left
h1=plot(Vlist,eta0/(1.44*8065.543887),'-','LineWidth',6, 'Color',[0,0,0]/255);hold on;
h2=plot(Vlist,etat/(1.44*8065.543887),'-','LineWidth',6, 'Color',[119,45,133]/255);hold on;
h3=plot(Vlist,eta2t/(1.44*8065.543887),'-','LineWidth',6, 'Color',[100,136,192]/255);hold on;
ylim([0 16])
ax.YColor = [0 0 0];
yticks([0 8 16]);
yticklabels( {'0','8','16'});
ylabel('P (units of $\omega_{21}$)','Interpreter','latex', 'FontSize', 24); 

yyaxis right
h4=plot(Vlist,100*etac,'LineWidth',6, 'Color',[238,90,57]/255);hold on;
ylim([-5 100])
yticks([0 50 95 100]);
yticklabels( {'0%','50%','','100%'});
ylabel('$\eta_{th}$','Interpreter','latex', 'FontSize', 28); 
ax.YColor = [238,90,57]/255;
ax.YColor = [0,0,0]/255;
ax.Box      = 'on';
ax.Layer    = 'top';
ax.LineWidth = 2;
ax.TickDir  = 'in';
set(ax, 'FontSize', 18);
% xlim(ax, [0 1]);
xticks([0 0.4 0.8 0.9 1 1.2]);
xticklabels({'0','0.4','0.8','0.9','1.0','1.2'});
xlabel('$V$ (units of $\omega_{21}$)', 'Interpreter', 'latex', 'FontSize', 24); 

lgd=legend([h1, h2, h3,h4], {'P($\Omega$=0)','P($\Omega$=$\Omega_{\rm T}$)','P($\Omega$=5$\Omega_{\rm T}$)','$\eta_{th}$'}, 'Interpreter', 'latex');
lgd.Position = [0.25, 0.58, 0.1, 0.1];  
lgd.Box = 'off';
lgd.FontSize = 22;

function [etac,eta0,etat,eta2t] = fun01(Vlist,nh,rh,rc,Gc)
nc=0;
Nc=0;
Eac=0.20*8065.543887;
Eab=1.82*8065.543887;
n=length(Vlist);
etac=zeros(1,n);
eta0=zeros(1,n);
etat=zeros(1,n);
eta2t=zeros(1,n);
Geff=(rh*(nh+1) + rc*(nc+1));
for i =1:n
    V=Vlist(i);
    G = Gc*exp((1.44*8065.543887)*(1-V)/208.5);
    Omeac=sqrt(0.5 * (G + rc*nc)*0.5 * (rh*(nh+1) + rc*(nc+1)));
    Omeac1=5*sqrt(0.5 * (G + rc*nc)*0.5 * (rh*(nh+1) + rc*(nc+1)));
    etac(1,i)=V*(1.44*8065.543887)/(1.82*8065.543887);
    eta0(1,i)=(G*Gc*nh*rc*rh*((417*log(Gc/G))/2 + (1.44*8065.543887)))/(G*Gc*rc + G*Gc*rh + 2.0*G*Gc*nh*rh + G*nh*rc*rh + Gc*nh*rc*rh);
    etat(1,i)=(G*Gc*nh*rh*((417*log(Gc/G))/2 + (1.44*8065.543887))*(G*rc + rc*rh + 4.0*Omeac^2 + rc^2 + nh*rc*rh))/(G^2*nh*rc*rh + 2.0*Gc*G^2*nh*rh + Gc*G^2*rc + Gc*G^2*rh + 4.0*G*Omeac^2*nh*rh + 4.0*Gc*G*Omeac^2 + G*nh^2*rc*rh^2 + 2.0*Gc*G*nh^2*rh^2 + G*nh*rc^2*rh + G*nh*rc*rh^2 + 4.0*Gc*G*nh*rc*rh + 3.0*Gc*G*nh*rh^2 + Gc*G*rc^2 + 2.0*Gc*G*rc*rh + Gc*G*rh^2 + 12.0*Gc*Omeac^2*nh*rh + 4.0*Gc*Omeac^2*rh + Gc*nh^2*rc*rh^2 + Gc*nh*rc^2*rh + Gc*nh*rc*rh^2);
    eta2t(1,i)=(G*Gc*nh*rh*((417*log(Gc/G))/2 + (1.44*8065.543887))*(G*rc + rc*rh + 4.0*Omeac1^2 + rc^2 + nh*rc*rh))/(G^2*nh*rc*rh + 2.0*Gc*G^2*nh*rh + Gc*G^2*rc + Gc*G^2*rh + 4.0*G*Omeac1^2*nh*rh + 4.0*Gc*G*Omeac1^2 + G*nh^2*rc*rh^2 + 2.0*Gc*G*nh^2*rh^2 + G*nh*rc^2*rh + G*nh*rc*rh^2 + 4.0*Gc*G*nh*rc*rh + 3.0*Gc*G*nh*rh^2 + Gc*G*rc^2 + 2.0*Gc*G*rc*rh + Gc*G*rh^2 + 12.0*Gc*Omeac1^2*nh*rh + 4.0*Gc*Omeac1^2*rh + Gc*nh^2*rc*rh^2 + Gc*nh*rc^2*rh + Gc*nh*rc*rh^2);
end
end
