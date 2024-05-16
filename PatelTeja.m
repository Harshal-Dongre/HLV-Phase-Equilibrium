function [fugacity,Z_vap] = PatelTeja(Pressure,Temperature,Pc,Tc,acentric)

% OUTPUT give [Vapour_Root Liquid_Root Fugacity] Function takes [Pressure Temperature] as input to the system

%% fugacity calculation by Patel Teja EOS
%%
% Patel Teja Parameters from acentric factor
F = 0.452413+1.30982*acentric-0.295937*acentric^2;
Zeta = 0.329032-0.076799*acentric+0.0211947*acentric^2;

R = 8.3144598;                             % Gas constant in [ MPa-cm3/(mol - K)]
Tr = Temperature/Tc;                                 % reduced Temp.

%  Omega Value Calculation - for parameter estimation of Patel Teja EOS
alpha = (1+F*(1-sqrt(Tr))).^2;             % Temp function
Omega_c = 1-3*Zeta;
Om_b = [1 (2-3*Zeta) (3*(Zeta^2)) -(Zeta^3)];   % Omega_b^3+(2-3*Zeta)*Omega_b^2+3*Zeta^2*Omega_b-Zeta^3 = 0;
Om_bs = roots(Om_b);
Om_bs_det = isreal(Om_bs);
if Om_bs_det == 1
    Omega_b = min(Om_bs);
else
    Omega_b = Om_bs(imag(Om_bs)==0);
end
Omega_a = 3*(Zeta^2)+3*(1-2*Zeta)*Omega_b+Omega_b^2+1-3*Zeta;
a = (Omega_a*alpha)*((R^2)*(Tc^2))/Pc;     %[MPa-m6/mol2]
b = Omega_b*R*Tc/Pc;                       %[m3/mol]
c = Omega_c*R*Tc/Pc;                       %[m3/mol]

% Compressibility Factor Estimation
A = (a*Pressure)/(R^2*(Temperature^2));                      %[unitless]
B = (b*Pressure)/(R*Temperature);                            %[unitless]
C = (c*Pressure)/(R*Temperature);                            %[unitless]

% maximum root for vapour phase fugacity calculations and minimum for liquid phase fugacity
Z_coeff = [1 (C-1) (A-2*B*C-B^2-B-C) (B*C+C-A)*B];   % Com_z = z^3+(C-1)*z^2+(A-2*B*C-B^2-B-C)*z+(B*C+C-A)*B;
Z_roots = roots(Z_coeff);
Z_vap = max(Z_roots(Z_roots >0 & imag(Z_roots) == 0));
Z_liq = min(Z_roots(Z_roots >0 & imag(Z_roots) == 0));

% Fugacity Coefficient Evaluation -Patel Teja Equation of state
N = (b*c+((b+c)^2)/2)^(-0.5);               %[cm3/mol]
M = (((b+c)/2)-N)*(Pressure/(R*Temperature));                %[unitless]
Q = (((b+c)/2)+N)*(Pressure/(R*Temperature));                %[unitless]

% molar volume of component calculations molar_volume_guest = Z_vap*R*Temperature/Pressure;

ln_osm_coeff_vap = Z_vap-1-log(Z_vap-B)+(a/(2*R*Temperature*N))*(log((Z_vap+M)/(Z_vap+Q)));       %[unitless] % fugacity coefficient vapour
ln_osm_coeff_liq = Z_liq-1-log(Z_liq-B)+(a/(2*R*Temperature*N))*(log((Z_liq+M)/(Z_liq+Q)));       %[unitless] % fugacity coefficient vapour
fugacity= exp([ln_osm_coeff_vap; ln_osm_coeff_liq])*Pressure;
end
