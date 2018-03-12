clc
clf
clear all;
close all;
 
%% Defining globals
global X T_0 M P_0 k0_1 Ea_1 R rho void_fraction U tube_CSA area_HX flow_area no_tubes Cv_steam m_steam Molar_mass catalyst_diameter  
%% 
%% Defining the variables
M_range = linspace(0,120000,500); % Defining the span of catalyst volume for which the solution is determined
T_0=(350+273.15);  % K
P_0=20*10^0 ; % bar
% k0_1= exp(7.4);% -
k0_1= 2.01*10^8;% -
%Activation energies for reaction 1 and reaction 2 (J/mol)
Ea_1= 114*1000; % J/mol
%Universal Gas Constant in J/(molK)
R=8.314;
rho=21.6; %?
void_fraction=0.5;
catalyst_diameter=0.0008;
%Saturated water is used as coolant at 5 bar.
Cv_steam=1588.16; 
m_steam=50; 

U=0.0;
%The diameter of the tubes are set such that 30 catalyst particles fit 
%perfectly without movement 
tube_OD=0.0130; % 3.50 inches
tube_ID=0.0127; % 3.15-inch, economical reasons easily available in market
tube_CSA=3.14*tube_ID^2/4; %cross sectional area of the tube
area_HX=4/tube_OD; %area avilable for heat exchange per unit volume of reactor (Fogler, Elements of Chemical Reaction Engineering)
no_tubes=350; %total number of tubes for the reactor
flow_area=no_tubes*tube_CSA; %total area for fluid flow
 
%% 
%% Defining the initial conditions---
m_H2_out=1.5; %kmol/s
Molar_mass=[2.02 28.01 18.02 44.01 16.04];
%Molar flow rate of cumene out in mol/s
F_H2_out=m_H2_out*1000;
% F_H2_in=1.283*1000;   
% F_CO_in=0.3016*1000;   
% F_H2O_in=1.0608*1000;   
% F_CO2_in=0.1635*1000;  
% F_CH4_in=0.00515*1000;  
F_H2_in=4146.300*1000/3600;    % mol/s
F_CO_in=1409.300*1000/3600;% mol/s
F_H2O_in=4090.700*1000/3600;   % mol/s
F_CO2_in=0*1000/3600;  % mol/s
F_CH4_in=32.300*1000/3600;  % mol/s
F_tot_in= F_H2_in+F_CO_in+F_H2O_in+F_CO2_in+F_CH4_in; % mol/s 
%The initial conditions can be arranged in a matrix z_0
z_0=[F_H2_in F_CO_in F_H2O_in F_CO2_in F_CH4_in F_tot_in T_0 P_0];
 
%% 
%% Now calculating the solution
[M,z] = ode45('HighTempODE', M_range, z_0);
%% 
%% The output matrix is z which has the following elements
F_H2=z(:,1);   
F_CO=z(:,2);   
F_H2O=z(:,3);  
F_CO2=z(:,4); 
F_CH4=z(:,5);  
F_tot=z(:,6);  
T_reactor=z(:,7);
P=z(:,8); 
%% Plotting the solution
%%
%Plotting the flowrates of the different compounds in the reactor
figure(1),...
    plot(M ,F_H2, '-r', M, F_CO,'-b', M, F_H2O, M, F_CO2,'-m', M, F_CH4),...
    title('PFR Flow Rates'), xlabel('Catalyst Weight, kg'),...
    ylabel('Flow rate, mol/s'),...
    legend('H2', 'CO', 'H2O', 'CO2', 'CH4');
 
%%
%Conversion of CO
X=(1-F_CO./F_CO_in);
 
%Plotting the conversion of propene against reactor volume
figure(2), plot(M,X, 'k--'),...
    title('Conversion of CO with reactor... volume'),...
    xlabel('Catalyst Weight, kg'), ylabel('Conversion, (-)');
 
%%
%Plotting the reactor pressure drop against reactor volume
figure(3), plot(M,P, 'k'), title('Pressure throughout the reactor'),...
    xlabel('Catalyst Weight, kg'), ylabel('Pressure, bar');
%%
%Plotting the coolant and reactor temp profile against reactor volume
figure(4), plot(M, T_reactor-273,'-b'),...
    title('Reactor Temperature Profile'),...
    xlabel('Catalyst Weight, kg'),ylabel('Temperature, degC'),...
    legend('T_reactor')

Hydrogen_Output= F_H2;
CO_Output= F_CO;
CO2_Output= F_CO2;
H2O_Output= F_H2O;


% %%
% %Finding the volume of the reactor 
Catalyst_weight=find(F_H2>1390);
G= find(M>0.1);
% %Finding volume of the catalyst required
Catalyst_volume=M(Catalyst_weight(1))/1250;
Reactor_volume=Catalyst_volume/void_fraction;
% 

formatSpec1='The volume of the required fixed bed tubular reactor is %1.3f m3\n';
fprintf(formatSpec1, Reactor_volume)
%  
formatSpec1='The volume of catalyst required is %1.3f m3\n';
fprintf(formatSpec1, Catalyst_volume)
%  
formatSpec1='The reactor conversion is about %1.2f \n';
fprintf(formatSpec1, X(Catalyst_weight(1)))
%  
formatSpec1='The reactor product temperature is %1.2f C\n';
fprintf(formatSpec1, T_reactor(Catalyst_weight(1))-273)
%  
delta_P=P_0/100000-P(Catalyst_weight(1));
formatSpec1='The tube-side pressure drop is %1.2f bar\n';
fprintf(formatSpec1, delta_P) 
% 
formatSpec1='The flow rate of Hydrogen is %1.2f kmol/h\n';
fprintf(formatSpec1, 3.6*F_H2(Catalyst_weight(1)))
% 
formatSpec1='The flow rate of Carbon Monoxide is %1.2f kmol/h\n';
fprintf(formatSpec1, 3.6*F_CO(Catalyst_weight(1)))
% 
formatSpec1='The flow rate of Carbon dioxide is %1.2f kmol/h\n';
fprintf(formatSpec1, 3.6*F_CO2(Catalyst_weight(1)))
% 
formatSpec1='The flow rate of Water is %1.2f kmol/h\n';
fprintf(formatSpec1, 3.6*F_H2O(Catalyst_weight(1)))
% 
formatSpec1='The flow rate of CH4 is %1.2f kmol/h\n';
fprintf(formatSpec1, 3.6*F_CH4(Catalyst_weight(1)))