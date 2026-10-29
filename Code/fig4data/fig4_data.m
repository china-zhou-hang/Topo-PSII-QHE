clear; clc;
%% ---------------- numerical controls ----------------
nOmega = 121;
nRc = 101;
nEps = 61;
nDeltaC = 35;
nPhiC = 25;
nOmegaMC = 81;
nSample = 700;

rng(20260723,'twister');

%% ---------------- basis and operators ----------------
e1 = [1;0;0;0];
e2 = [0;1;0;0];
e3 = [0;0;1;0];
e4 = [0;0;0;1];

P1 = e1*e1';
P2 = e2*e2';
P3 = e3*e3';
P4 = e4*e4';

S14 = e1*e4';  S14d = S14';
S34 = e3*e4';  S34d = S34';
S12 = e1*e2';  S12d = S12';
S23 = e2*e3';  S23d = S23';

%% ---------------- reference parameters ----------------
eV2cm = 8065.543887;
kB = 0.6950;
Tc = 300;
Th = 6000;
p0.E1 = 0;
p0.E2 = 0.18*eV2cm;
p0.E3 = 1.62*eV2cm;
p0.E4 = 1.82*eV2cm;
p0.omega43 = p0.E4 - p0.E3;
p0.Eab = p0.E3 - p0.E2;

p0.nh = 10000;
p0.nc = 0;
p0.Nc = 0;

p0.rh = 2.23e-7*eV2cm;
p0.rc = 2.23e-3*eV2cm;
p0.Gc = 0.0248*eV2cm;
p0.etaC = 1-Tc/Th;

p0.V = 0.95;
p0.G = p0.Gc*exp(p0.Eab*(1-p0.V)/(kB*Tc));

p0.Omega = 0;
p0.Delta = 0;
p0.gamma_phi = 0;

p0.P1 = P1; p0.P2 = P2; p0.P3 = P3; p0.P4 = P4;
p0.S14 = S14; p0.S14d = S14d;
p0.S34 = S34; p0.S34d = S34d;
p0.S12 = S12; p0.S12d = S12d;
p0.S23 = S23; p0.S23d = S23d;

GammaD0 = p0.rh*(p0.nh+1) + p0.rc*(p0.nc+1);
GammaA0 = p0.G + p0.rc*p0.nc;
OmegaT0 = 0.5*sqrt(GammaD0*GammaA0);

fprintf('Reference parameters:\n');
fprintf('  G = %.6g\n',p0.G);
fprintf('  Gamma_D = %.6g\n',GammaD0);
fprintf('  Gamma_A = %.6g\n',GammaA0);
fprintf('  Omega_T/omega_43 = %.6f\n',OmegaT0/p0.omega43);

%% ========================================================================
%  Part I. Panels (a) and (b): retained unchanged
%  ========================================================================
OmegaList = linspace(0,0.6,nOmega)*p0.omega43;
rcFactorList = linspace(0.1,1.5,nRc);

etaMap = zeros(nRc,nOmega);
powerMap = zeros(nRc,nOmega);
JMap = zeros(nRc,nOmega);
cohMap = zeros(nRc,nOmega);
nuMap = nan(nRc,nOmega);

fprintf('\nComputing (Omega,r_c) map for panels (a) and (b)...\n');
for iO = 1:nOmega
    for iR = 1:nRc
        p = p0;
        p.Omega = OmegaList(iO);
        p.rc = p0.rc*rcFactorList(iR);

        obs = qhe_steady_state_rwa(p);

        etaMap(iR,iO) = obs.eta;
        powerMap(iR,iO) = obs.power;
        JMap(iR,iO) = obs.Jcoh;
        cohMap(iR,iO) = obs.coherence;
        nuMap(iR,iO) = effective_response_winding(p);
    end

    if mod(iO,max(1,round(nOmega/10))) == 0 || iO == nOmega
        fprintf('  map progress: %d/%d\n',iO,nOmega);
    end
end

% Transverse logarithmic sensitivity to r_c.
sensitivityMap = zeros(size(etaMap));
logRc = log(rcFactorList(:));
for iO = 1:nOmega
    y = log(max(etaMap(:,iO),1.0e-14));
    sensitivityMap(:,iO) = abs(gradient(y,logRc));
end

mask0 = (nuMap == 0);
mask1 = (nuMap == 1);
meanS0 = mean(sensitivityMap(mask0),'omitnan');
meanS1 = mean(sensitivityMap(mask1),'omitnan');

fprintf('\nMean transverse logarithmic sensitivity:\n');
fprintf('  nu_eff = 0: %.6g\n',meanS0);
fprintf('  nu_eff = 1: %.6g\n',meanS1);
fprintf('  ratio S_0/S_1 = %.6g\n',meanS0/max(meanS1,eps));

%% ========================================================================

OmegaC = 1.20*OmegaT0;
DeltaScale = 0.60*p0.omega43;
GammaPhiScale = 1.50*GammaD0;
epsilonList = linspace(0,1.50,nEps);

pC0 = p0;
pC0.Omega = OmegaC;
obsC0 = qhe_steady_state_rwa(pC0);
etaC0 = obsC0.eta;
JC0 = max(abs(obsC0.Jcoh),eps);

etaWorst = ones(1,nEps);
JWorst = ones(1,nEps);
marginWorst = nan(1,nEps);
nuFraction = nan(1,nEps);

fprintf('\nComputing fixed-drive worst-case robustness for panel (c)...\n');
for iE = 1:nEps
    epsNow = epsilonList(iE);

    DeltaListC = linspace(-epsNow*DeltaScale, ...
                          epsNow*DeltaScale,nDeltaC);
    PhiListC = linspace(0,epsNow*GammaPhiScale,nPhiC);

    etaLocal = inf;
    JLocal = inf;
    marginLocal = inf;
    nNu1 = 0;
    nTot = 0;

    for iP = 1:numel(PhiListC)
        for iD = 1:numel(DeltaListC)
            p = p0;
            p.Omega = OmegaC;
            p.Delta = DeltaListC(iD);
            p.gamma_phi = PhiListC(iP);

            obs = qhe_steady_state_rwa(p);
            [nuNow,marginNow] = effective_response_winding(p);

            etaLocal = min(etaLocal,obs.eta/etaC0);
            JLocal = min(JLocal,abs(obs.Jcoh)/JC0);
            marginLocal = min(marginLocal,marginNow/OmegaT0);

            nNu1 = nNu1 + double(nuNow == 1);
            nTot = nTot + 1;
        end
    end

    etaWorst(iE) = etaLocal;
    JWorst(iE) = JLocal;
    marginWorst(iE) = marginLocal;
    nuFraction(iE) = nNu1/nTot;

    if mod(iE,max(1,round(nEps/10))) == 0 || iE == nEps
        fprintf('  worst-case progress: %d/%d\n',iE,nEps);
    end
end

% Locate the first closing of the response point gap for the entire
% uncertainty set: min(m_top)=0.
epsilonT = NaN;
idClose = find(marginWorst <= 0,1,'first');
if ~isempty(idClose)
    if idClose == 1
        epsilonT = epsilonList(1);
    else
        e1c = epsilonList(idClose-1);
        e2c = epsilonList(idClose);
        m1c = marginWorst(idClose-1);
        m2c = marginWorst(idClose);
        epsilonT = e1c + (0-m1c)*(e2c-e1c)/(m2c-m1c);
    end
end

fprintf('\nPanel (c):\n');
fprintf('  fixed Omega/omega_43 = %.6f\n',OmegaC/p0.omega43);
fprintf('  point-gap-closing uncertainty radius epsilon_T = %.6f\n',epsilonT);

%% ========================================================================


mc.sigmaOmega = 0.06;                       % 6% rms drive fluctuation
mc.sigmaRc = 0.10;                          % 10% log-normal r_c disorder
mc.sigmaG = 0.10;                           % 10% log-normal extraction disorder
mc.deltaMax = 0.25*p0.omega43;              % |Delta| <= 0.25 omega_43
mc.gammaPhiMax = 1.20*GammaD0;              % 0 <= gamma_phi <= 1.2 Gamma_D

OmegaMC = linspace(0,0.6,nOmegaMC)*p0.omega43;

etaNominal = zeros(1,nOmegaMC);
ratioMean = zeros(1,nOmegaMC);
ratioMedian = zeros(1,nOmegaMC);
ratioP05 = zeros(1,nOmegaMC);
ratioP95 = zeros(1,nOmegaMC);
ratioCV = zeros(1,nOmegaMC);
probNu1 = zeros(1,nOmegaMC);
marginMean = zeros(1,nOmegaMC);

fprintf('\nComputing combined-disorder robustness for panel (d)...\n');
for iO = 1:nOmegaMC
    pNom = p0;
    pNom.Omega = OmegaMC(iO);
    obsNom = qhe_steady_state_rwa(pNom);
    etaNominal(iO) = obsNom.eta;

    ratioSamples = zeros(nSample,1);
    nuSamples = nan(nSample,1);
    marginSamples = nan(nSample,1);

    for k = 1:nSample
        p = p0;

        p.Omega = max(0,OmegaMC(iO)*(1+mc.sigmaOmega*randn));
        p.Delta = mc.deltaMax*(2*rand-1);
        p.rc = p0.rc*exp(mc.sigmaRc*randn);
        p.G = p0.G*exp(mc.sigmaG*randn);
        p.gamma_phi = mc.gammaPhiMax*rand;

        obs = qhe_steady_state_rwa(p);
        ratioSamples(k) = obs.eta/max(etaNominal(iO),eps);

        [nuSamples(k),marginSamples(k)] = effective_response_winding(p);
    end

    ratioMean(iO) = mean(ratioSamples);
    ratioMedian(iO) = median(ratioSamples);
    ratioP05(iO) = percentile_no_toolbox(ratioSamples,5);
    ratioP95(iO) = percentile_no_toolbox(ratioSamples,95);
    ratioCV(iO) = std(ratioSamples)/max(mean(ratioSamples),eps);
    probNu1(iO) = mean(nuSamples == 1);
    marginMean(iO) = mean(marginSamples,'omitnan')/OmegaT0;

    if mod(iO,max(1,round(nOmegaMC/10))) == 0 || iO == nOmegaMC
        fprintf('  Monte-Carlo progress: %d/%d\n',iO,nOmegaMC);
    end
end

%% ========================================================================
xOmega = OmegaList/p0.omega43;
xOmegaMC = OmegaMC/p0.omega43;

save('fig4_data.mat','p0','OmegaT0','xOmega','xOmegaMC','rcFactorList','etaMap','sensitivityMap','nuMap','epsilonList','etaWorst','JWorst','marginWorst','epsilonT','ratioP05','ratioP95','ratioMedian','probNu1');

%% ========================================================================
%  Local functions
%  ========================================================================

function obs = qhe_steady_state_rwa(p)
% Solve the complete four-level RWA Lindblad steady state by vectorization.

    N = 4;
    I = eye(N);

    H = -p.Delta*p.P4 + p.Omega*(p.S34+p.S34d);

    L = -1i*(kron(I,H)-kron(H.',I));

    L = L + p.rh*((p.nh+1)*dissipator_super(p.S14) ...
                           + p.nh*dissipator_super(p.S14d));

    L = L + p.rc*((p.nc+1)*dissipator_super(p.S34) ...
                           + p.nc*dissipator_super(p.S34d));

    L = L + p.Gc*((p.Nc+1)*dissipator_super(p.S12) ...
                           + p.Nc*dissipator_super(p.S12d));

    L = L + p.G*dissipator_super(p.S23);

    % Pure-dephasing convention:
    % (gamma_phi/2)D[P4-P3] gives d rho_43/dt=-gamma_phi*rho_43.
    if p.gamma_phi > 0
        Z = p.P4-p.P3;
        L = L + 0.5*p.gamma_phi*dissipator_super(Z);
    end

    traceRow = zeros(1,N^2);
    for j = 1:N
        traceRow(j+(j-1)*N) = 1;
    end

    A = L;
    b = zeros(N^2,1);
    A(end,:) = traceRow;
    b(end) = 1;

    rhoVec = A\b;
    rho = reshape(rhoVec,N,N);

    % Remove numerical non-Hermiticity and normalize.
    rho = 0.5*(rho+rho');
    rho = rho/real(trace(rho));

    rho33 = real(rho(3,3));       % acceptor population
    rho43 = rho(4,3);             % donor-acceptor coherence

    obs.rho = rho;
    obs.rho33 = rho33;
    obs.coherence = 2*abs(rho43);
    obs.Jcoh = 2*p.Omega*imag(rho43);
    obs.power = p.G*rho33*(p.V*p.Eab);
    obs.eta = obs.power/(p.etaC*p.nh*p.rh*(p.E4-p.E1));
end

function DS = dissipator_super(A)
% Column-major vectorization:
% vec(A*rho*A^dagger)=(A^* kron A)vec(rho).

    N = size(A,1);
    I = eye(N);
    AdA = A'*A;

    DS = kron(conj(A),A) ...
       - 0.5*kron(I,AdA) ...
       - 0.5*kron(AdA.',I);
end

function [nu,margin,znum,zden] = effective_response_winding(p)

    GammaD = p.rh*(p.nh+1)+p.rc*(p.nc+1);
    GammaA = p.G+p.rc*p.nc;

    a = 0.5*GammaD+0.25*p.gamma_phi;
    b = 0.5*GammaA+0.25*p.gamma_phi;

    c1 = a+b-1i*p.Delta;
    c0 = b*(a-1i*p.Delta);

    znum = roots([1,c1,c0-p.Omega^2]);
    zden = roots([1,c1,c0+p.Omega^2]);

    scale = max([1;abs(znum(:));abs(zden(:))]);
    tol = 1.0e-9*scale;

    if any(abs(real(znum)) <= tol) || any(abs(real(zden)) <= tol)
        nu = NaN;
    else
        nu = sum(real(znum)>tol)-sum(real(zden)>tol);
        nu = round(real(nu));
    end

    margin = max(real(znum));
end

function q = percentile_no_toolbox(x,pct)
% Percentile with linear interpolation, without Statistics Toolbox.

    x = sort(x(:));
    n = numel(x);

    if n == 0
        q = NaN;
        return;
    elseif n == 1
        q = x;
        return;
    end

    u = 1+(n-1)*(pct/100);
    i1 = floor(u);
    i2 = ceil(u);

    if i1 == i2
        q = x(i1);
    else
        q = x(i1)+(u-i1)*(x(i2)-x(i1));
    end
end
