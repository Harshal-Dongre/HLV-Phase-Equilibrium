function vdWP_difference = vdWP_minimizer(T_exp,P_exp,guest_conc,salt_wt,Tc,Pc,Vc,acentric,Henry_Parameters,UNIQUAC_identifier,...
    SRK_identifier,Kihara_potential,Guest_atomic_properties,salt_identifier,UNIQUAC_data,hydrate_structure)

%%
% % Predicted pressures - empty storage matrix
vdWP_difference = deal(zeros(numel(T_exp),7));
P_Pred_tmp = deal(zeros(numel(T_exp),1));

structure = hydrate_structure;              % hydrate structure allotment

for u = 1:numel(T_exp)
    T = T_exp(u);

    % structure properties and hydrate surface properties

    if structure == 1
        cg_pr_wtr = [1/23 3/23];                    % [sI - number of cages of type per water molecute]
        if T > 273.15
            del_hw0 = -4858.9;                      % [J/mol]
            VP_constants = [4.1539 -5500.9332 7.6537 -16.1277e-3];          % [WATER] => molar volume [m3/mol] Klauda Sandler 2000
        elseif T < 273.15
            del_hw0 = 1151;                         % [J/mol]
            VP_constants = [4.6056 -5501.1243 2.9446 -8.1431e-3];           % [ICE]
        end
        del_vw0 = 4.6;                              % [cm3/mol]
        del_mu_bulk_0 = 1263.60;                    % [J/mol]
    elseif structure == 2
        cg_pr_wtr = [2/17 1/17];                    % [sII - number of cages of type per water molecute]
        if T > 273.15
            VP_constants = [4.1539 -5500.9332 7.6537 -16.1277e-3];          % [WATER]
            del_hw0 = -5202.2;                      % [J/mol]
        elseif T < 273.15
            del_hw0 = 808;                          % [J/mol]
            VP_constants = [4.6056 -5501.1243 2.9446 -8.1431e-3];           % [ICE] => KS 2000 [m3/mol]
        end
        del_vw0 = 5.0;                              % [cm3/mol]
        del_mu_bulk_0 = 882.80;                     % [J/mol]
    end
    
    P_sat_wtr = exp(VP_constants(1)*log(T)+VP_constants(2)/T+VP_constants(3)+VP_constants(4)*T)*1e-6 ;   % [MPa]    

    % determining the composition of the system
    if numel(Tc) > 1                  % represents mixed guest formers
        guest_fraction = guest_conc(u,:);
    elseif numel(Tc) == 1               % represents single and pure hydrate former
        guest_fraction = [];
        SRK_identifier = [];
        SRK_interactions = [];
    end
    
    if numel(salt_wt) >= 1 & salt_wt ~= 0
        salt_wt_percent = salt_wt(u,:);
        % Pitzer activity coefficient model
        if numel(salt_wt_percent) == 1
            ln_gamma_salt = Pitzer_single(salt_wt_percent,salt_identifier);
        else
            ln_gamma_salt = Pitzer_mixed(salt_wt_percent_in,salt_identifier);
        end
        
        % Setschnow salting out factor and the number of ions per hydrate lattice
        [salting_out_factor,ion_per_lattice,component_mol_frac] = salt_conversion ...
            (salt_wt_percent, salt_identifier, Kihara_potential(3),ideal_water_molecules);
        
    else
        salt_wt_percent = [];
        ln_gamma_salt = 0;
        salting_out_factor = 1;
        ion_per_lattice = 0;
        component_mol_frac = 0;
    end
    
    % temperature based global estimates
    Henry_Constant = (exp(-(Henry_Parameters(:,1)+Henry_Parameters(:,2)./T+ ...
        Henry_Parameters(:,3).*log(T)+Henry_Parameters(:,4).*T))*0.101325)';  % [MPa]
    
    % molar volume of guests in water at infinite dilution, molar_volume_inf => [Zhou and Battino 2001 10.1021/je000215o]
    guest_mlr_vlm_inf = (10.74 + 0.2683*Vc*1000);     % [cm3/mol]

    % Kihara potential estimate
    Langmuir_const = Kihara_mix(T,structure,Kihara_potential);

    P = P_exp(u)*[0.5 1.5];
    iter = 0;
    
    for i = 1:400
        [err_P1,Z1_vap,vapor_fugacity1] = vdWP_estimate(T,P(1),guest_fraction,Tc,Pc,acentric,cg_pr_wtr,del_hw0,del_vw0,del_mu_bulk_0,...
            P_sat_wtr,Henry_Constant,guest_mlr_vlm_inf,Langmuir_const,...
            UNIQUAC_identifier,UNIQUAC_data,SRK_identifier,SRK_interactions);
        [err_P2,Z2_vap,vapor_fugacity2] = vdWP_estimate(T,P(2),guest_fraction,Tc,Pc,acentric,cg_pr_wtr,del_hw0,del_vw0,del_mu_bulk_0,...
            P_sat_wtr,Henry_Constant,guest_mlr_vlm_inf,Langmuir_const,...
            UNIQUAC_identifier,UNIQUAC_data,SRK_identifier,SRK_interactions);
        err = [err_P1 err_P2];
        iter = iter+1;

        if prod(err)<0
            P_itr = mean(P);
            [err_itr,Zitr_vap,vapor_fugacityitr] = vdWP_estimate(T,P_itr,guest_fraction,Tc,Pc,acentric,cg_pr_wtr,del_hw0,del_vw0,del_mu_bulk_0,...
                P_sat_wtr,Henry_Constant,guest_mlr_vlm_inf,Langmuir_const,...
                UNIQUAC_identifier,UNIQUAC_data,SRK_identifier,SRK_interactions);
            if err(1)*err_itr < 0
                P(2) = P_itr;
            else
                P(1) = P_itr;
            end
            P_Pred_tmp(u) = P_itr;
            Z_Pred_tmp(u) = Zitr_vap;
            vapor_fugacity_tmp(u) = vapor_fugacityitr;

            % termination criterion
            if abs(diff(err)) < 1e-7
                break
            end
        else
            if P > 1
                P(1) = (P(1))^(1-iter/100);
                P(2) = (P(2))^(1+iter/100);
            else
                P(1) = (P(1))^(1+iter/100);
                P(2) = (P(2))^(1-iter/100);
            end
        end
        % garbage value
        if iter > 300  || P(2) > 5*P_exp(u) || P(1) < 1e-5
            break
        end
    end
    Y = [T P_exp(u) guest_conc(u) P_Pred_tmp(u) Z_Pred_tmp(u) vapor_fugacity_tmp(u) structure];    % dummy output matrix
    % disp(Y)
    disp([T P_exp(u) P_Pred_tmp(u)])
    vdWP_difference(u,:) = Y;
end

end
