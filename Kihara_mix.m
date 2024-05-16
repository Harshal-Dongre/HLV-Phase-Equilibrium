function C_MPa_out = Kihara_mix(Temperature,structure,potential_parameters)

%%
%{
potential parameters are in the order of CORE_RADIUS | SIGMA | EPSILON per row and different cloumns indicate different guests
%}

% Constants
k_B = 1.38064852e-23;                               % [J/K] => Boltzmann constant

if structure == 1
    Cavity_Radius = [3.95e-10 4.33e-10];        % [m]                             ===>> Sloan Book
    Z_LJp = [20 24];                            % coordination number             ===>> Sloan Book - sI type hydrate
elseif structure == 2
    Cavity_Radius = [3.91e-10 4.73e-10];        % [m]                             ===>> Sloan Book
    Z_LJp = [20 28];                            % coordination number             ===>> Sloan Book - sII type hydrate
end

C_MPa_out = deal(zeros(numel(potential_parameters(1,:)),2));   % zeros allocations for storage
[C_Pa C_MPa] = deal(zeros(1,2));


for i = 1:numel(potential_parameters(1,:))
    Core_radius = potential_parameters(1,i)*1e-10;
    sigma = potential_parameters(2,i)*1e-10;
    epsilon = potential_parameters(3,i);
    % delN is a distance function for cavities
    for j =1:2
        del = @(r,N) ((1/N)*((1-(r./Cavity_Radius(j))-(Core_radius/Cavity_Radius(j))).^(-N)-(1+(r./Cavity_Radius(j))-(Core_radius/Cavity_Radius(j))).^(-N)));      % ==> Original function
        % potential function for the Langmuir constants
        W = @(r) (2.*Z_LJp(j)*epsilon*(((sigma^12./(Cavity_Radius(j).^11.*r)).*(del(r,10)+Core_radius.*del(r,11)./Cavity_Radius(j)))-((sigma^6./(Cavity_Radius(j).^5.*r)).*(del(r,4)+Core_radius.*del(r,5)./Cavity_Radius(j)))));
        C_trm = @(r) (exp(-W(r)/(Temperature))).*(r.^2);
        C_Pa(j) = 4*pi.*(integral(C_trm,0,(Cavity_Radius(j)-Core_radius)))./(k_B.*Temperature); %===>> Unit ==>> [1/Pa]  % Langmuir constants evaluated
        C_MPa(j) =  C_Pa(j)*1e6;                                    % ==>> Inverse relationship ==> multiply by (10^6) instead of (10^-6)
        C_MPa_out(i,j) = C_MPa(j);
    end
end

end
