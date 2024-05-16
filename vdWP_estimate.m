function [error,Z_vap,vapor_fugacity] = vdWP_estimate(T,P,guest_fraction,Tc,Pc,acentric,cg_pr_wtr,del_hw0,del_vw0,del_mu_bulk_0,...
            P_sat_wtr,Henry_Constant,guest_mlr_vlm_inf,Langmuir_const,...
            UNIQUAC_identifier,UNIQUAC_data,SRK_identifier,SRK_interactions)      
        
R = 8.3144598;                      % [(cm3.MPa)/(mol.K)]
        
%% chemical potential estimate of hydrate lattice

% Fugacity estimate from the Patel Teja EoS model
if numel(Tc) == 1
    [fugacity_PT,Z_vap] = PatelTeja(P,T,Pc,Tc,acentric);
elseif numel(Tc) > 1
    [fugacity_PT,Z_vap] = Patel_Teja_Mixed(T,P,guest_fraction,Tc,Pc,acentric,SRK_identifier,SRK_interactions);          % [vapor and liquid respectively];
end
vapor_fugacity = fugacity_PT(1,:);
liquid_fugacity = fugacity_PT(2,:);

% fractional occupancy
theta = ((Langmuir_const'.*vapor_fugacity))'./(1+(vapor_fugacity*Langmuir_const));

% chemical potential hypothetical empty hydrate lattice
del_mu_hyd = -sum(cg_pr_wtr.*log(1-sum(theta,1)));     % [cm3.MPa/mol]

%% chemical potential estimate of the equilibrium water phase

% mole fraction calculations
x_guest = liquid_fugacity./(Henry_Constant.*exp(guest_mlr_vlm_inf.*(P-P_sat_wtr)/(R*T)));     % mole fraction of guest in pure water
x_water = 1-sum(x_guest);                   % mole fraction of water

% activity coefficient calculations
ln_gamma = UNIQUAC(UNIQUAC_data,UNIQUAC_identifier,x_guest,T);
gamma = exp(ln_gamma(1));                   % activity coefficient of water

% Chemical potential of hydrate phase

del_hw = del_hw0 + integral (@(T) (-38.12 + 0.141 * (T-273.15)),273.15,T);  % [J/mol]
del_vw = @(P) (del_vw0 + 6.695e-15*P*1e6);
del_mu_bulk = del_mu_bulk_0/(R*273.15) ...
              - integral (@(T) del_hw./(R*T.^2),273.15,T) ...
              + integral (@(P) del_vw(P)/(R*T),0,P) - log(gamma*x_water);

% potential difference between the chemical potential estimates of the equilibrium phases
error = del_mu_bulk - del_mu_hyd;
end
