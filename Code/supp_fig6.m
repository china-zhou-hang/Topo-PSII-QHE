clear all;clc;
load('fig2data/cmap_RdYlBu.mat');
cmap = cmap_RdYlBu;
set(groot,'DefaultFigureColormap',cmap);

%% -------------------- parameters --------------------
p.eV2cm = 8065.543887;

p.E0  = 0;
p.E1  = 0.18*p.eV2cm;
p.E2  = 1.62*p.eV2cm;
p.E3  = 1.82*p.eV2cm;

p.E21 = 1.44*p.eV2cm;
p.E32 = 0.20*p.eV2cm;

p.nh      = 1.0e4;
p.gamma_h = 2.23e-7*p.eV2cm;
p.gamma_c = 2.23e-3*p.eV2cm;
p.Gamma_c = 0.0248*p.eV2cm;

p.kB   = 0.6950;
p.Tc   = 300;
p.kBTc = p.kB*p.Tc;

%% -------------------- voltage scan --------------------
Vlist = linspace(0.8,1.05,200);

Omega_T     = zeros(size(Vlist));
Omega_Tstar = zeros(size(Vlist));
Omega_TL    = zeros(size(Vlist));
Omega_c     = zeros(size(Vlist));

for k = 1:numel(Vlist)
    S = transport_scales(Vlist(k),p);

    Omega_T(k)     = S.Omega_T;
    Omega_Tstar(k) = S.Omega_Tstar;
    Omega_TL(k)    = S.Omega_TL;

    % Independently extract Omega_c from the maximum of the exact
    % steady-state donor-acceptor coherence. The optimization variable
    % x = Omega/Omega_TL avoids numerical difficulties when the scale
    % changes by many orders of magnitude across the voltage scan.
    Omega_c(k) = coherence_peak_scaled(S,p);
end

%% -------------------- diagnostics --------------------
relErr_TL = abs(Omega_TL-Omega_c)./Omega_c;
relErr_Ts = abs(Omega_Tstar-Omega_c)./Omega_c;
save('supp_fig6_data.mat','Vlist','p','Omega_T','Omega_Tstar','Omega_TL','Omega_c','relErr_TL','relErr_Ts');
%% 
close all
color1=[218,243,241]/255;
color2=[253,235,208]/255;
figure('Position', [-1000,200,700,400]);
ax1= subplot(1,1,1);
patch(ax1, [0.9 1.1 1.1 0.9], [-1 -1 10 10],color1, 'FaceAlpha',1, 'EdgeColor','none');hold on;
patch(ax1, [1.05 1.0 1.0 1.05], [-1 -1 10 10],color2, 'FaceAlpha',1, 'EdgeColor','none');hold on;
h1=plot(Vlist,Omega_T/p.E32,'-','LineWidth',4, 'Color',[0,0,0]/255);hold on;
h2=plot(Vlist,Omega_Tstar/p.E32,'-','LineWidth',4, 'Color',[238,90,57]/255);hold on;
h3=plot(Vlist,Omega_TL/p.E32,'-.','LineWidth',2, 'Color',[100,136,192]/255);hold on;
% plot(Vlist,Omega_c/p.E32,':','LineWidth',4, 'Color',[238,90,57]/255);hold on;
idxMark = 1:10:numel(Vlist);
h4=semilogy(Vlist(idxMark),Omega_c(idxMark)/p.E32,'o','LineStyle','none','Color',[0.12 0.38 0.68],'MarkerFaceColor','w','MarkerSize',5.8,'LineWidth',1.2);
plot(0.95, 0.1, 'p', 'MarkerSize', 16, 'MarkerEdgeColor', [238,90,57]/255, 'MarkerFaceColor', [255,255,255]/255,'LineWidth', 1.8);
ax1.Box      = 'on';
ax1.Layer    = 'top';
ax1.LineWidth = 2;
ax1.TickDir  = 'in';
set(ax1, 'FontSize', 20);
xlim([0.8 1.05])
xticks([0.8 0.9 0.95 1.0 1.05]);
xticklabels({'0.8','0.9','0.95','1.0','1.05'});
ylim([-0.5 9])
yticks([0 3 6 9]);
yticklabels({'0','3.0','6.0','9.0'});
xlabel('$V$ (units of $\omega_{21}$)', 'Interpreter', 'latex', 'FontSize', 24); 
ylabel('critical drive$/\omega_{32}$','Interpreter','latex', 'FontSize', 24); 
% text(0.45,0.15,'(0.95,1.27)', 'Units', 'normalized', 'FontSize', 24,Color=[0,0,0])
lgd=legend([h1, h2,h3,h4], {'$\Omega_T$','$\Omega_T^{*}$','$\Omega_{T,\mathcal{L}}$','$\Omega_c$'}, 'Interpreter', 'latex');
lgd.Position = [0.55, 0.57, 0.1, 0.1];  
lgd.Box = 'off';
lgd.FontSize = 22;
lgd.NumColumns = 2;
%% =====================================================
% Local functions
% ======================================================

function S = transport_scales(V,p)
% Calculate the four-level transport rates and critical scales at fixed
% normalized voltage V/omega_21.

% Extraction rate:
% Gamma = Gamma_c exp[-(V_phys-omega_21)/(k_B T_c)]
% with V_phys = V*omega_21.
S.Gamma = p.Gamma_c*exp(p.E21*(1-V)/p.kBTc);

S.Rh = p.gamma_h*(p.nh+1);
S.Wh = p.gamma_h*p.nh;

% Unified notation used throughout Supplementary Note 3.
S.Gamma_eff = S.Rh+p.gamma_c;
S.kappa = (S.Gamma+S.Gamma_eff)/2;

S.n1 = S.Gamma+2*p.Gamma_c+S.Rh+2*S.Wh;

S.n0 = S.Gamma*p.Gamma_c...
     + S.Gamma*S.Wh...
     + p.Gamma_c*S.Rh...
     + 2*p.Gamma_c*S.Wh;

S.d2 = S.Gamma+p.Gamma_c+S.Rh+S.Wh+p.gamma_c;

S.d1 = S.Gamma*p.Gamma_c...
     + (S.Gamma+p.Gamma_c)*(S.Rh+S.Wh+p.gamma_c)...
     + S.Wh*p.gamma_c;

S.d0 = S.Gamma*p.Gamma_c*(S.Rh+S.Wh+p.gamma_c)...
     + S.Wh*p.gamma_c*(S.Gamma+p.Gamma_c);

% Bare no-jump threshold.
S.Omega_T = 0.5*sqrt(S.Gamma*S.Gamma_eff);

% Pade-type feedback-renormalized approximation.
S.alpha = ((3*p.Gamma_c+S.Wh)*S.Gamma...
         + S.Wh*(p.Gamma_c+S.Wh))...
         /(2*(S.Gamma+S.Wh)*(p.Gamma_c+S.Wh));

S.Gamma_eff_star = S.alpha*S.Gamma_eff;
S.Omega_Tstar = 0.5*sqrt(S.Gamma*S.Gamma_eff_star);

% Exact complete-Liouvillian transport-response threshold.
S.Omega_TL = sqrt(S.kappa*S.d0/(2*S.n0));
end

function C = exact_coherence(Omega,S,p)
% Exact steady-state l1-norm donor-acceptor coherence in the resonant
% RWA transport sector.

numerator = 2*p.Gamma_c*S.Wh...
          * abs(S.Gamma-p.gamma_c).*Omega;

denominator = S.kappa*S.d0...
            + 2*S.n0.*Omega.^2;

C = numerator./denominator;
end

function Omega_c = coherence_peak_scaled(S,p)
% Numerically determine Omega_c from max C_l1.
%
% The dimensionless optimization variable x=Omega/Omega_TL keeps the
% search stable over the full voltage interval.

objective = @(x) -exact_coherence(x*S.Omega_TL,S,p);

opts = optimset(...
    'TolX',1.0e-12,...
    'Display','off');

xPeak = fminbnd(objective,0,5,opts);
Omega_c = xPeak*S.Omega_TL;
end
