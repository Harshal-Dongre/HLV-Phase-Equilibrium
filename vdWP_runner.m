function vdWP_difference = vdWP_runner(data,constants,number_components,UNIQUAC_data)
%% Caller function and data split

% experimental data
T_exp = data(:,1);
P_exp = data(:,2);

% the number of components are matricized as number of guests and number of salt species. Thus, if only one guest
%is present it has to be assumed to be 100% guest concentration (pure)
if number_components(1)==1
    guest_conc = ones(size(data,1),1);                                                  % represents pure species
    salt_wt = data(:,(3:(2+number_components(2))));
else
    guest_conc = data(:,(3:(2+number_components(1))));        % guest concentration in vapor phase [y-mole frac]
    salt_wt = data(:,(2+number_components(1)+1:(2+sum(number_components))));
end

% critical constants of hydrate formers
Tc = constants(1,:);
Pc = constants(2,:);
Vc = constants(3,:);
acentric = constants(4,:);

% Henry Parameters temperature dependent
Henry_Parameters = (constants(5:8,1:numel(Tc)))';

% UNIQUAC identifier
UNIQUAC_identifier = constants(9,:);

% Gas phase SRK interactions
SRK_identifier = constants(10,:);

% Kihara potential parameters
Kihara_potential = constants(11:13,1:numel(Tc));

% Non covalent interaction constants
Guest_NCI = constants(15:17,1:numel(Tc));

% Salt constants  - salting out estimates
salt_identifier = constants(18);

% Hydrate structure identifier
hydrate_structure = constants(19);

% SRK interaction parameters for mixed guest systems not included yet.
vdWP_difference =  vdWP_minimizer(T_exp,P_exp,guest_conc,salt_wt ,Tc,Pc,Vc,acentric,Henry_Parameters,UNIQUAC_identifier,...
    SRK_identifier,Kihara_potential,Guest_NCI,salt_identifier,UNIQUAC_data,hydrate_structure);

end
