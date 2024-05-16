function ln_gamma = UNIQUAC(UNIQUAC_data,identifier,guest_fraction,Temperature)

%{
***** DONT FORGET TO LOAD THE UNIQUAC_data.mat ***** INPUT NEEDED - UNIQUAC interaction matrix - ensure same order UNIQUAC pure parameters - ensure same order identifier - ensure same order
guest_fraction Temperature
%}

%%              ***** ----- ACTIVITY OF WATER IS THE FIRST VALUE OF THE MATRIX ----- *****
%%      IDENTIFIER   ||      COMPONENT
%{
        1           ||      Water
        2           ||      Carbon dioxide [CO2]
        3           ||      Oxygen [O2]
        4           ||      Nitrogen [N2]
        5           ||      Methane [CH4]
        6           ||      Ethane [C2H6]
        7           ||      Propane [C3H8]
%}
%%
%% INTERACTION PARAMETER IDENTIFICATION

UNIQUAC_interaction_matrix = UNIQUAC_data{1,1};
UNIQUAC_pure_parameters = UNIQUAC_data{1,2};

% component water incorporated in the system
identifier = [1 identifier];                     % water identifier included
identifier = sort(identifier);                   % setting in order

% vectors of identifiers for combination
combinations_2_make = {identifier,identifier};
interaction_elements = sortrows(combvec(combinations_2_make{:}).');

% linear index of elements mapped to the indexes of the interaction database
interaction_index = sub2ind(size(UNIQUAC_interaction_matrix),interaction_elements(:,2),...
    interaction_elements(:,1));

interaction_isolation = UNIQUAC_interaction_matrix(interaction_index);   % isolate required values from the original matrix

% interactions needed for the model
interaction = reshape(interaction_isolation,numel(identifier),numel(identifier));

%% UNIQUAC activity coefficient model

% gas law constant
R = 8.3144598;

r = UNIQUAC_pure_parameters(identifier,2)';
q = UNIQUAC_pure_parameters(identifier,1)';

% concentration matrix
x = [1-sum(guest_fraction) guest_fraction];

% UNIQUAC parameters phi = x.*r/sum(x.*r);
theta = x.*q/sum(x.*q);

% temperature dependent interaction
tau = exp(-interaction/(R*Temperature));

% residual and combinatorial function evaluation
J = r./sum(x.*r);
L = q./sum(x.*q);
s = theta*tau;

% residuals contribution ln_gamma_residual = q(1)*(1-log(s(1))-(theta./s)*tau(1,:)')       % first row only - water
ln_gamma_residual = q.*(1-log(s)-(theta./s)*tau');                 % full matrix

% combinatorial contribution
ln_gamma_combination = 1 - J + log (J) - 5*q.* (1- ((J./L) + log (J./L)));

% activity coefficient of guests
ln_gamma = ln_gamma_residual+ln_gamma_combination;
end
