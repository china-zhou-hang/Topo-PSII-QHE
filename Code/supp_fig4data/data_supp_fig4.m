clear all;close all;clc;
V=0.95;
eV2cm=8065.543887;
G=(0.0248*eV2cm)*exp((1.44*eV2cm)*(1-V)/(0.6950*300)); 
rh = 2.23e-7*eV2cm;
rc = 2.23e-3*eV2cm;
nh=10000;
nc = 0;
Gc=0.0248*eV2cm;
a3 = 0.5 * (G + rc*nc);
a4 = 0.5 * (rh*(nh+1) + rc*(nc+1));
Ome_T = sqrt(a3*a4);
%% 
t0=0.01/5305;
tf=4000/5305;
Omega1 = 0;
Omega2 = 0.5*Ome_T;
Omega3 = 1*Ome_T;
Omega4 = 5*Ome_T;
Omega5 = 10*Ome_T;

tlist=logspace(log10(t0),log10(tf),1000);
% tlist=linspace(t0,tf,1000);
% tlist=logspace(-4,0,3000);
% [t_1, p44_1, p33_1,p43_1] = fun01([t0 tf], G, Omega1,V);
% [t_2, p44_2, p33_2,p43_2] = fun01([t0 tf], G, Omega2,V);
% [t_3, p44_3, p33_3,p43_3] = fun01([t0 tf], G, Omega3,V);
[t_1, p44_1, p33_1, p22_1, p11_1,p43_1,p4p3_1,Jc_1] = fun01(tlist, G, Omega1,nh,rh,rc,Gc);
[t_2, p44_2, p33_2, p22_2, p11_2,p43_2,p4p3_2,Jc_2] = fun01(tlist, G, Omega2,nh,rh,rc,Gc);
[t_3, p44_3, p33_3, p22_3, p11_3,p43_3,p4p3_3,Jc_3] = fun01(tlist, G, Omega3,nh,rh,rc,Gc);
[t_4, p44_4, p33_4, p22_4, p11_4,p43_4,p4p3_4,Jc_4] = fun01(tlist, G, Omega4,nh,rh,rc,Gc);
[t_5, p44_5, p33_5, p22_5, p11_5,p43_5,p4p3_5,Jc_5] = fun01(tlist, G, Omega5,nh,rh,rc,Gc);
tlist=5305*tlist;

save('supp_fig4_data.mat','tlist','p44_*','p33_*','p22_*','p11_*','p43_*','p4p3_*');
%% 
figure('Position', [-1700,50,1600,1200]);
subplot(2,2,1)
plot(t_1*5305, p44_1,'r-', 'LineWidth',1.5); hold on;
plot(t_2*5305, p44_2,'g-', 'LineWidth',1.5); hold on;
plot(t_3*5305, p44_3,'b-', 'LineWidth',1.5); hold on;
plot(t_3*5305, p44_4,'y-', 'LineWidth',1.5); hold on;
plot(t_3*5305, p44_5,'k-', 'LineWidth',1.5); hold on;
set(gca, 'XScale', 'log');
% xlim([0,200])

subplot(2,2,2)
plot(t_1*5305, p33_1,'r-', 'LineWidth',1.5); hold on;
plot(t_2*5305, p33_2,'g-', 'LineWidth',1.5); hold on;
plot(t_3*5305, p33_3,'b-', 'LineWidth',1.5); hold on;
plot(t_3*5305, p33_4,'y-', 'LineWidth',1.5); hold on;
plot(t_3*5305, p33_5,'k-', 'LineWidth',1.5); hold on;
set(gca, 'XScale', 'log');
% xlim([0,200])

subplot(2,2,3)

plot(t_1*5305, p4p3_1,'r-', 'LineWidth',1.5); hold on;
plot(t_2*5305, p4p3_2,'g-', 'LineWidth',1.5); hold on;
plot(t_3*5305, p4p3_3,'b-', 'LineWidth',1.5); hold on;
plot(t_3*5305, p4p3_4,'y-', 'LineWidth',1.5); hold on;
plot(t_3*5305, p4p3_5,'k-', 'LineWidth',1.5); hold on;
set(gca, 'XScale', 'log');

subplot(2,2,4)
plot(t_1*5305, Jc_1,'r-', 'LineWidth',1.5); hold on;
plot(t_2*5305, Jc_2,'g-', 'LineWidth',1.5); hold on;
plot(t_3*5305, Jc_3,'b-', 'LineWidth',1.5); hold on;
plot(t_3*5305, Jc_4,'y-', 'LineWidth',1.5); hold on;
plot(t_3*5305, Jc_5,'k-', 'LineWidth',1.5); hold on;
set(gca, 'XScale', 'log');

%% 
function [t, p44, p33,p22,p11,p43,p4p3,Jc] = fun01(tspan, G, Omeac,nh,rh,rc,Gc)
    e1=[1;0;0;0]; e2=[0;1;0;0]; e3=[0;0;1;0]; e4=[0;0;0;1];
    P1=e1*e1'; P2=e2*e2'; P3=e3*e3'; P4=e4*e4';
    S14=e1*e4'; S14d=S14';
    S34=e3*e4'; S34d=S34';
    S12=e1*e2'; S12d=S12';
    S23=e2*e3'; S23d=S23';
    E41=(1.82*8065.543887)/(1.44*8065.543887);
    nc=0; Nc=0;
    phi=0;
    p=struct('P1',P1,'P2',P2,'P3',P3,'P4',P4, ...
             'S14',S14,'S14d',S14d,'S34',S34,'S34d',S34d, ...
             'S12',S12,'S12d',S12d,'S23',S23,'S23d',S23d, ...
             'rh',rh,'rc',rc,'Gc',Gc,'G',G,'nh',nh,'nc',nc,'Nc',Nc, ...
             'Omeac',Omeac,'phi',phi);
    y0 = P1(:);
    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t,y] = ode15s(@(tt,yy) rhs_rwa(tt,yy,p), tspan, y0, opts);
    nT = numel(t);
    p44 = zeros(nT,1);
    p4p3     = zeros(nT,1);
    p43     = zeros(nT,1);
    p33     = zeros(nT,1);
    p22     = zeros(nT,1);
    p11     = zeros(nT,1);
    Jc     = zeros(nT,1);
    for k=1:nT
        R = reshape(y(k,:),4,4);
        p44(k) = real(R(4,4));
        p43(k) = 2*abs(R(4,3));
        p33(k) = real(R(3,3));
        p22(k) = real(R(2,2));
        p11(k) = real(R(1,1));
        p4p3(k) = real(R(4,4))-real(R(3,3));
        Jc(k)=2*Omeac*abs(R(4,3));
    end

end
function dy = rhs_rwa(~,y,p)
    rho = reshape(y,4,4);
    % HRWA = p.Omega*(exp(1i*p.phi)*(p.S34d + p.S12d) + exp(-1i*p.phi)*(p.S34 + p.S12)) ...
    %    - p.Delta*p.P4 - (p.Delta + 144)*p.P2;
        HRWA = p.Omeac*(exp(1i*p.phi)*(p.S34d) + exp(-1i*p.phi)*(p.S34));
    drho = -1i*(HRWA*rho - rho*HRWA);
    drho = drho + dissipators(rho,p);
    dy = drho(:);
end

function L = D(rho,A)
    L = A*rho*A' - 0.5*(A'*A*rho + rho*A'*A);
end

function drho = dissipators(rho,p)
    drho = zeros(4,4);
    drho = drho + p.rh*((p.nh+1)*D(rho,p.S14) + p.nh*D(rho,p.S14d));
    drho = drho + p.rc*((p.nc+1)*D(rho,p.S34) + p.nc*D(rho,p.S34d));
    drho = drho + p.Gc*((p.Nc+1)*D(rho,p.S12) + p.Nc*D(rho,p.S12d));
    drho = drho + (p.G/2)*(2*p.S23*rho*p.S23d - p.S23d*p.S23*rho - rho*p.S23d*p.S23);
end
