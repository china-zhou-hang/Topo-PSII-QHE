function value = function_fig2(mode,varargin)
% Functions used by data_fig2.m. All formulae retain their original forms.
switch mode
    case 'rhs_norwa'
        value = rhs_norwa(varargin{:});
    case 'rhs_rwa'
        value = rhs_rwa(varargin{:});
    case 'rwa_eta'
        value = rwa_eta(varargin{:});
    case 'rwa_hessian'
        value = rwa_hessian(varargin{:});
    case 'rwa_rho32'
        value = rwa_rho32(varargin{:});
    otherwise
        error('Unknown function_fig2 mode: %s',mode);
end
end

%% ---------- RWA analytical expressions ----------
function value = rwa_eta(Omeac,G,nh,rh,rc,Gc,E21,E30,kBTc)
value = (20*G*Gc*(E21 + kBTc*log(Gc/G))*(G*rc + rc*rh + 4.0*Omeac^2 + rc^2 + nh*rc*rh))/(19*E30*(G^2*nh*rc*rh + 2.0*Gc*G^2*nh*rh + Gc*G^2*rc + Gc*G^2*rh + 4.0*G*Omeac^2*nh*rh + 4.0*Gc*G*Omeac^2 + G*nh^2*rc*rh^2 + 2.0*Gc*G*nh^2*rh^2 + G*nh*rc^2*rh + G*nh*rc*rh^2 + 4.0*Gc*G*nh*rc*rh + 3.0*Gc*G*nh*rh^2 + Gc*G*rc^2 + 2.0*Gc*G*rc*rh + Gc*G*rh^2 + 12.0*Gc*Omeac^2*nh*rh + 4.0*Gc*Omeac^2*rh + Gc*nh^2*rc*rh^2 + Gc*nh*rc^2*rh + Gc*nh*rc*rh^2));
end

function value = rwa_hessian(Omeac,G,nh,rh,rc,Gc,E21,E30,kBTc)
value = (160*G*Gc*(E21 + kBTc*log(Gc/G)))/(19*E30*(G^2*nh*rc*rh + 2*Gc*G^2*nh*rh + Gc*G^2*rc + Gc*G^2*rh + 4*G*Omeac^2*nh*rh + 4*Gc*G*Omeac^2 + G*nh^2*rc*rh^2 + 2*Gc*G*nh^2*rh^2 + G*nh*rc^2*rh + G*nh*rc*rh^2 + 4*Gc*G*nh*rc*rh + 3*Gc*G*nh*rh^2 + Gc*G*rc^2 + 2*Gc*G*rc*rh + Gc*G*rh^2 + 12*Gc*Omeac^2*nh*rh + 4*Gc*Omeac^2*rh + Gc*nh^2*rc*rh^2 + Gc*nh*rc^2*rh + Gc*nh*rc*rh^2)) + (40*G*Gc*(E21 + kBTc*log(Gc/G))*(8*G*Gc*Omeac + 8*Gc*Omeac*rh + 8*G*Omeac*nh*rh + 24*Gc*Omeac*nh*rh)^2*(G*rc + rc*rh + 4*Omeac^2 + rc^2 + nh*rc*rh))/(19*E30*(G^2*nh*rc*rh + 2*Gc*G^2*nh*rh + Gc*G^2*rc + Gc*G^2*rh + 4*G*Omeac^2*nh*rh + 4*Gc*G*Omeac^2 + G*nh^2*rc*rh^2 + 2*Gc*G*nh^2*rh^2 + G*nh*rc^2*rh + G*nh*rc*rh^2 + 4*Gc*G*nh*rc*rh + 3*Gc*G*nh*rh^2 + Gc*G*rc^2 + 2*Gc*G*rc*rh + Gc*G*rh^2 + 12*Gc*Omeac^2*nh*rh + 4*Gc*Omeac^2*rh + Gc*nh^2*rc*rh^2 + Gc*nh*rc^2*rh + Gc*nh*rc*rh^2)^3) - (20*G*Gc*(E21 + kBTc*log(Gc/G))*(8*Gc*rh + 8*G*Gc + 8*G*nh*rh + 24*Gc*nh*rh)*(G*rc + rc*rh + 4*Omeac^2 + rc^2 + nh*rc*rh))/(19*E30*(G^2*nh*rc*rh + 2*Gc*G^2*nh*rh + Gc*G^2*rc + Gc*G^2*rh + 4*G*Omeac^2*nh*rh + 4*Gc*G*Omeac^2 + G*nh^2*rc*rh^2 + 2*Gc*G*nh^2*rh^2 + G*nh*rc^2*rh + G*nh*rc*rh^2 + 4*Gc*G*nh*rc*rh + 3*Gc*G*nh*rh^2 + Gc*G*rc^2 + 2*Gc*G*rc*rh + Gc*G*rh^2 + 12*Gc*Omeac^2*nh*rh + 4*Gc*Omeac^2*rh + Gc*nh^2*rc*rh^2 + Gc*nh*rc^2*rh + Gc*nh*rc*rh^2)^2) - (320*G*Gc*Omeac*(E21 + kBTc*log(Gc/G))*(8*G*Gc*Omeac + 8*Gc*Omeac*rh + 8*G*Omeac*nh*rh + 24*Gc*Omeac*nh*rh))/(19*E30*(G^2*nh*rc*rh + 2*Gc*G^2*nh*rh + Gc*G^2*rc + Gc*G^2*rh + 4*G*Omeac^2*nh*rh + 4*Gc*G*Omeac^2 + G*nh^2*rc*rh^2 + 2*Gc*G*nh^2*rh^2 + G*nh*rc^2*rh + G*nh*rc*rh^2 + 4*Gc*G*nh*rc*rh + 3*Gc*G*nh*rh^2 + Gc*G*rc^2 + 2*Gc*G*rc*rh + Gc*G*rh^2 + 12*Gc*Omeac^2*nh*rh + 4*Gc*Omeac^2*rh + Gc*nh^2*rc*rh^2 + Gc*nh*rc^2*rh + Gc*nh*rc*rh^2)^2);
end

function value = rwa_rho32(Omeac,G,nh,rh,rc,Gc)
value = (2*Gc*Omeac*nh*rh*(G*1.0i - rc*1.0i))/(G^2*nh*rc*rh + 2.0*Gc*G^2*nh*rh + Gc*G^2*rc + Gc*G^2*rh + 4.0*G*Omeac^2*nh*rh + 4.0*Gc*G*Omeac^2 + G*nh^2*rc*rh^2 + 2.0*Gc*G*nh^2*rh^2 + G*nh*rc^2*rh + G*nh*rc*rh^2 + 4.0*Gc*G*nh*rc*rh + 3.0*Gc*G*nh*rh^2 + Gc*G*rc^2 + 2.0*Gc*G*rc*rh + Gc*G*rh^2 + 12.0*Gc*Omeac^2*nh*rh + 4.0*Gc*Omeac^2*rh + Gc*nh^2*rc*rh^2 + Gc*nh*rc^2*rh + Gc*nh*rc*rh^2);
end

%% ---------- Lindblad master equations ----------
function dy = rhs_rwa(t,y,p)
rho = reshape(y,4,4);
H_RWA = -p.Delta*p.P3 + p.Omega*(p.S23 + p.S23d);
drho = -1i*(H_RWA*rho - rho*H_RWA);
drho = drho + dissipators(rho,p);
dy = drho(:);
end

function dy = rhs_norwa(t,y,p)
rho = reshape(y,4,4);
wdrive = p.omega32 + p.Delta;
Hdrive = 2*p.Omega*cos(wdrive*t + p.phi)*(p.S23 + p.S23d);
Hlab = p.H0 + Hdrive;
drho = -1i*(Hlab*rho - rho*Hlab);
drho = drho + dissipators(rho,p);
dy = drho(:);
end

function L = D(rho,A)
L = A*rho*A' - 0.5*(A'*A*rho + rho*A'*A);
end

function drho = dissipators(rho,p)
drho = zeros(4,4);
drho = drho + p.rh*((p.nh+1)*D(rho,p.S03) + p.nh*D(rho,p.S03d));
drho = drho + p.rc*((p.nc+1)*D(rho,p.S23) + p.nc*D(rho,p.S23d));
drho = drho + p.Gc*((p.Nc+1)*D(rho,p.S01) + p.Nc*D(rho,p.S01d));
drho = drho + (p.G/2)*(2*p.S12*rho*p.S12d- p.S12d*p.S12*rho-rho*p.S12d*p.S12 );
end
