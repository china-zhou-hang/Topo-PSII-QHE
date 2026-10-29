% data_fig2.m
% Calculate and save all numerical and analytical data for Figure 2.
clear; clc;
tic

%% ---------- basis & operators ----------
e0=[1;0;0;0]; e1=[0;1;0;0]; e2=[0;0;1;0]; e3=[0;0;0;1];
P0=e0*e0'; P1=e1*e1'; P2=e2*e2'; P3=e3*e3';
S03=e0*e3'; S03d=S03';
S23=e2*e3'; S23d=S23';
S01=e0*e1'; S01d=S01';
S12=e1*e2'; S12d=S12';

%% ---------- static parameters ----------
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
phi=0;

p.P0=P0; p.P1=P1; p.P2=P2; p.P3=P3;
p.S03=S03; p.S03d=S03d; p.S23=S23; p.S23d=S23d;
p.S01=S01; p.S01d=S01d; p.S12=S12; p.S12d=S12d;
p.rh=rh; p.rc=rc; p.Gc=Gc; p.nh=nh; p.nc=nc; p.Nc=Nc;
p.phi=phi; p.H0=H0; p.omega32=E32;
V=0.95;
p.Delta=0;
p.G=Gc*exp(E21*(1-V)/kBTc);
y0=P0(:);

%% ---------- scan ranges ----------
Omelist = linspace(0,0.6,100)*E32;
Omelist1 = linspace(0,0.6,10)*E32;
Glist = logspace(log10(1e1),log10(1e6),100);
nO=numel(Omelist);
nG=numel(Glist);

%% ---------- parallel pool: 50 workers ----------
pool = gcp('nocreate'); if ~isempty(pool), delete(pool); end; parpool('threads',50);

%% ---------- Fig. 2(a:i): Non-RWA 2D scan ----------
wdrive = p.omega32+p.Delta;
T=2*pi/abs(wdrive);
nPeriods_relax=300;
nPeriods_avg=20;
PtsPerPeriod=400;
PtsPerPeriod_avg=40;
Tend_relax=(nPeriods_relax+nPeriods_avg)*T;
tspan_relax=[0,(nPeriods_relax-1)*T,nPeriods_relax*T+(0:nPeriods_avg*PtsPerPeriod_avg)*(T/PtsPerPeriod_avg)];
opts=odeset('RelTol',1e-12,'AbsTol',1e-12,'MaxStep',T/PtsPerPeriod);

[IO,JG]=ndgrid(1:nO,1:nG);
IO=IO(:);
JG=JG(:);
Ngrid=numel(IO);
Exact_eta_lin=zeros(Ngrid,1);

parfor n=1:Ngrid
    iO=IO(n);
    jG=JG(n);
    p_loc=p;
    p_loc.Omega=Omelist(iO);
    p_loc.G=Glist(jG);
    [tNR,yNR]=ode15s(@(t,y) function_fig2('rhs_norwa',t,y,p_loc),tspan_relax,y0,opts);

    selN=tNR >= nPeriods_relax*T;
    tAvg=tNR(selN);y_last_N=yNR(selN,:);
    Nn=size(y_last_N,1);
    rho22=zeros(Nn,1);
    for k=1:Nn
        R=reshape(y_last_N(k,:),4,4);
        rho22(k)=real(R(3,3));
    end
    Vloc=(kBTc*log(Gc/p_loc.G)+E21)/E21;
    Exact_eta_lin(n)=p_loc.G*(trapz(tAvg,rho22)/(tAvg(end)-tAvg(1)))*Vloc/(etaC*p_loc.nh*p_loc.rh*(E3-E0))*E21;
end
Exact_eta_2D=reshape(Exact_eta_lin,[nO,nG]).';

%% ---------- Fig. 2(a:ii): Non-RWA 1D scan (10 points) ----------
Tend=(nPeriods_relax+nPeriods_avg)*T;
Nt=(nPeriods_relax+nPeriods_avg)*PtsPerPeriod;
tspan=linspace(0,Tend,Nt);
Exact_eta_1D=zeros(1,numel(Omelist1));

parfor iO=1:numel(Omelist1)
    p_loc=p;
    p_loc.Omega=Omelist1(iO);
    [tNR,yNR]=ode15s(@(t,y) function_fig2('rhs_norwa',t,y,p_loc),tspan,y0,opts);
    selN=tNR >= (Tend-nPeriods_avg*T);
    tAvg=tNR(selN);y_last=yNR(selN,:);
    NN=size(y_last,1);
    rho22=zeros(NN,1);
    for k=1:NN
        R=reshape(y_last(k,:),4,4);
        rho22(k)=real(R(3,3));
    end
    Vloc=(kBTc*log(Gc/p_loc.G)+E21)/E21;
    Exact_eta_1D(iO)=p_loc.G*(trapz(tAvg,rho22)/(tAvg(end)-tAvg(1)))*Vloc/(etaC*p_loc.nh*p_loc.rh*(E3-E0))*E21;
end

%% ---------- Fig. 2(a:ii): RWA efficiency and Hessian (100 points) ----------
G=p.G;
RWA_eta_1D=zeros(1,nO);
RWA_detadO2_1D=zeros(1,nO);
for iO=1:nO
    Omeac=Omelist(iO);
    RWA_eta_1D(iO)=function_fig2('rwa_eta',Omeac,G,nh,rh,rc,Gc,E21,E30,kBTc);
    RWA_detadO2_1D(iO)=function_fig2('rwa_hessian',Omeac,G,nh,rh,rc,Gc,E21,E30,kBTc);
end

%% ---------- Fig. 2(c:i): coherence map ----------
Cl1_2D=zeros(nG,nO);
for j=1:nG
    Gj=Glist(j);
    for iO=1:nO
        Omeac=Omelist(iO);
        rho32=function_fig2('rwa_rho32',Omeac,Gj,nh,rh,rc,Gc);
        Cl1_2D(j,iO)=2*abs(rho32);
    end
end

%% ---------- Fig. 2(c:ii): coherence and coherent flux ----------
RWA_Cl1_1D=zeros(1,nO);
RWA_Jc_1D=zeros(1,nO);
for iO=1:nO
    Omeac=Omelist(iO);
    rho32=function_fig2('rwa_rho32',Omeac,G,nh,rh,rc,Gc);
    RWA_Cl1_1D(iO)=2*abs(rho32);
    RWA_Jc_1D(iO)=2*Omeac*imag(rho32);
end

%% ---------- save data ----------
if ~exist('fig2data','dir'), mkdir('fig2data'); end
save(fullfile('fig2data','fig2_data.mat'), ...
    'Exact_eta_2D','Exact_eta_1D','RWA_eta_1D','RWA_detadO2_1D', ...
    'Cl1_2D','RWA_Cl1_1D','RWA_Jc_1D');
toc
