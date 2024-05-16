%% Code for the estimation of the Hydrate-Liquid-Vapor phase equilibrium line.

%{
If you need help in understanding the application and use of this code, kindly go through the ReadMe file. 
Additionally, it also presents link to a brief discussion about this HLVE. 

Here, 
* Vapor phase ---> Patel Teja Equation of state
* Liquid phase ---> UNIQUAC activity coefficient model (guests)
* Hydrate phase ---> van der Waals and Platteeuw model 

A total of five single guest hydrates are presented as test case scenario. 
These can be initiated through the use of "select" variable.

Select Value        Guest
1                   Methane
2                   Ethane
3                   Propane
4                   Carbon dioxide
5                   R152a
%}

clc
clear
format shortG

%% Sub-model datasets
load UNIQUAC_data
load EoS_interactions

%% Excel datasheets
datasets_ranges = readtable('data_range.csv');

select = 5;
selected_range = datasets_ranges(select,:)

datasets_stringed = string(table2cell(selected_range));

% Experimental hydrate datasets
data = xlsread('HLVE_Experimental_data.xlsx','experiment',datasets_stringed(4));
constants = xlsread('HLVE_Experimental_data.xlsx','experiment',datasets_stringed(7));
number_components = str2double(datasets_stringed(5:6));         % [guests salts]

% Function executions
vdWP_difference = vdWP_runner(data,constants,number_components,UNIQUAC_data);

% AARD(P) calculations
Pexp = data(:,2);
P_Pred = vdWP_difference(:,4);

figure(1)
scatter(data(:,1),data(:,2),'r*','DisplayName','Experiment')
hold on
plot(data(:,1),P_Pred(:,1),'bO--','DisplayName','Model')
legend('Experiment', 'Model', 'Location','northwest')

AADP = (1/size(P_Pred,1))*sum(abs(Pexp-P_Pred)./Pexp,1)*100 % AADP
