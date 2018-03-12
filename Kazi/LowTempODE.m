%reactor_cumene.m
%
%Functional file of differential equations for sizing cumene reactor
%% 
% File written by Asko Kukk, Kazi Naveed, Stephanos, Intasar,Siyao
%% 
% Last modified 
 
%% Variable dictionary
%% Reaction Kinetics
%k0_i      - global -   Reaction constants for reaction 1 & 2(l/min h)
%Ea_i      - global -   The activation energy for reaction 1 & 2 (J/mol)
%r1        - variable - Reaction rate for main reaction
%r2        - variable - Reaction rate for side reaction
%% Reaction Enthalpy
% delta_H1  - global -   enthalpy change for reaction 1 (J/mol)
% delta_H2  - global -   enthalpy change for reaction 2 (J/mol)
%% 
% V_range   - Reactor_volumeor -   Reactor_volumeor defining volume range of catalyst required (m^3)
%% 
% Ratm      - global -   ideal gas constant in kcal/(molK)
% R         - global -   ideal gas constant in (J/mol K)
%% Catalyst Properties 
% rho               - global -   ideal gas constant in (J/mol K)
% void_fraction     - global -   ideal gas constant in (J/mol K)
% catalyst_diameter - global -   catalyst particle diameter (m)
%% Reactor tube dimensions
%tube_ID    - variable - Reactor tubes inner diameter (m)
%tube_OD    - variable - Reactor tubes outer diameter (m)
%tube_csa   - variable - Cross-sectional area of the tube (m2)
%area_HX    - variable - Heat transfer area per reactor volume (1/m)
%no_tubes   - variable - Number of reactor tubes
%flow_area  - variable - Total area for flow
%% Coolant Properties
%T_coolant_0 - global - Initial temperature of coolant (K)
%T_coolant   - global - Temperature of coolant in reactor (K)
%Cv_steam    - global - Heat Capacity of coolant at 5 bar (J/kgK)
%m_steam     - global - Mass flowrate of coolant
%% Heat Transfer Coeffcient
%U           - global - Overall heat transfer coefficient (W/m2K)
%% Reactor Initial Conditions
%T0            - variable -  Reactor feed temperature (K)
%P0            - variable -  Reactor inlet pressure (bar)
%F_cum_in      - variable -  Cumene feed flow rate (mol/s)
%F_ben_in      - variable -  Benzene feed flow rate (mol/s)
%F_propene_in  - variable -  Propene feed flow rate (mol/s)
%F_DIPB_in     - variable -  DIPB feed flow rate (mol/s)
%F_propane_in  - variable -  Propanefeed flow rate (mol/s)
%F_tot_in      - variable -  Total molar feed flow rate (mol/s)
%z_0           - variable -  Vector containing initial conditions of reactor

%% Reactor Outlet Conditions
%m_cum_out   - variable -  Required mass flow rate of cumene (kg/s)
%F_cum_out   - variable -  Required molar flow rate of cumene (mol/s)
%
%F_cum       - variable -  Cumene flow rate in reactor (mol/s)
%F_ben       - variable -  Benzene flow rate in reactor(mol/s)
%F_propene   - variable -  Propene flow rate in reactor(mol/s)
%F_DIPB      - variable -  DIPB flow rate in reactor(mol/s)
%F_propane   - variable -  Propane flow rate in reactor(mol/s)
%F_tot       - variable -  Total molar feed flow rate in reactor(mol/s)
%T_reactor   - variable -  Tube side reactor temperature (K)
%T_coolant   - variable -  Temperature of Coolant (K)
%P           - variable -  Tube-side pressure (bar)
%z           - variable -  Vector containing reactor conditions
%% 
%t           - variable -    Time / independent variable (s)
%X           - variable -    Conversion of Propene (h)
%V           - variable -    Volume of the reactor (m^3)
%zf          - variable -    Compressibility factor
%mu          - variable -    Viscosity
%Cp          - variable -    Heat capacity
 
function dzdM=ODE(V, z)
 
%Global variables from Cumene_reactor.m
global T_0 P_0 r1 k0_1 Ea_1 R rho void_fraction U tube_CSA area_HX delta_H1 delta_H2 Cv_steam m_steam Molar_mass catalyst_diameter flow_area no_tubes V_molar
F_H2=z(1);   %mol/s
F_CO=z(2); %mol/s  
F_H2O=z(3);  %mol/s
F_CO2=z(4); %mol/s
F_CH4=z(5);  %mol/s
F_tot=z(6);  %mol/s
T_reactor=z(7); %K
P=z(8); %bar
%% Defining the gas properties changing with temperature
%Specific heat capacities of the gases with respect to temperature
tt=T_reactor/1000;
Cp_H2=33.066178-11.36*tt+11.43*tt^2-2.77*tt^3-0.159/tt^2;
Cp_CO=25.56759+6.096*tt+4.05*tt^2-2.67*tt^3+0.131/tt^2;
Cp_H2O=30.09+6.833*tt+6.793*tt^2-2.53*tt^3 +0.08/tt^2;
Cp_CO2=24.997+55.19*tt-33.69*tt^2+7.95*tt^3-0.14/tt^2;
Cp_CH4=-0.703+108.48*tt-42.52*tt^2+5.86*tt^3+0.68/tt^2;
%Viscosities of the gases with respect to temperature
mu_H2=1.201*10^-5;
mu_CO=2.49*10^-5;
mu_H2O=0.000012;
mu_CO2=2.269*10^-5;
mu_CH4=1.578*10^-5;

%Weighted average viscosity of the gases
mu=F_H2/F_tot*mu_H2+F_CO/F_tot*mu_CO+F_H2O/F_tot*mu_H2O+F_CO2/F_tot*mu_CO2+F_CH4/F_tot*mu_CH4;

%Molar flow rate multiplied by molar volume gives volumetric flow rate
%l/s=(mol/s * l/mol)

%% Defining the rate constants k1 and k2 and reaction enthalpies
%k1= k0_1 *exp(-1*Ea_1/(R*T_reactor));
%k1= exp(-29364/(1.987*T_reactor)+40.32/1.987); %mol/gcat.min
k1= 1.85*10^-5*exp(12.88-1855.5/T_reactor);
dCp1std = 3.22087; % J/mol K basis h2o
delta_H1_std=-4.1*10^4; %J/mol
dCp1 = Cp_H2 + Cp_CO2 - Cp_CO - Cp_H2O; % J/mol K basis h2o
delta_H1 = delta_H1_std + (dCp1+dCp1std)/2*(T_reactor-(298.15)); %J/mol Heat of reaction (reference 25degC)

%% 
% The compressibility of gases is a big issue in the reactor. Therefore
% an average compressibility factor is assumed of a value of 0.9. This
% value is the average taken from the Unisim file. The compressibility
% factor is hence used to calculate the molar volume and volumetric
% flowrates
zf=0.9;

% Molar volume from the equation of state
V_molar=zf*R*T_reactor/P*1000; 
 
% Gas volumetric flowrate from molar volume
volumetric_flowrate=z(1:6).*V_molar;

% Calculating mass flow rates
mass_flowrate=z(1:5).*Molar_mass'; 

% Calculating total mass flow rate
m_tot=sum(mass_flowrate);

% Calculating superficial velocity
u=(volumetric_flowrate(6)/no_tubes)/(tube_CSA*1000);

%% Calculating gas density using density=mass/volume
rho_gas=(m_tot)/volumetric_flowrate(6); 

Kp=exp((4577.8/T_reactor)-4.33);

y_H2 = (F_H2 / F_tot) ; % Pa
y_CO = (F_CO / F_tot) ;
y_H2O = (F_H2O / F_tot) ;
y_CO2 = (F_CO2 / F_tot) ; 
y_CH4 = (F_CH4 / F_tot) ; 
k11= 454*exp(12.88-1855.5/T_reactor)/(0.0283*3600);
AF= 0.86 + 0.14*P;
rhobee= 1442;
r1= k11*AF*(y_H2O*y_CO-y_CO2*y_H2/Kp)/(rhobee*379);
%% Differential equations solving
% Cumene flow rate with respect to reactor volume 
dzdM(1)=r1;
% Benzene flow rate with respect to reactor volume 
dzdM(2)=-r1;
% Propene flow rate with respect to reactor volume 
dzdM(3)=-r1;
% DIPB flow rate with respect to reactor volume 
dzdM(4)=r1;
% Propane flow rate with respect to reactor volume 
dzdM(5)=0;
% Total flow rate with respect to reactor volume
dzdM(6)=dzdM(1)+dzdM(2)+dzdM(3)+dzdM(4)+dzdM(5);

% Pressure with respect to reactor volume in Pascals 
% The Ergun equation was used for the following equation:
dzdM(8)=-1*0*((150*mu*u*(1-void_fraction)^2)/(catalyst_diameter^2*void_fraction^3)+(1.75*rho_gas*u^2*(1-void_fraction))/(catalyst_diameter*void_fraction^3))/flow_area;

% Reactor temperature with respect to reactor volume
% Reactor Temperature with co-current coolant for PBR was used: 
% dzdV(7)=(U*area_HX*(T_coolant-T_reactor)+((-r1)*(delta_H1)))/(F_H2*Cp_H2+F_CO*Cp_CO+F_H2O*Cp_H2O+F_CO2*Cp_CO2+F_CH4*Cp_CH4);
dzdM(7)=(-1*r1)*(delta_H1)/(F_H2*Cp_H2+F_CO*Cp_CO+F_H2O*Cp_H2O+F_CO2*Cp_CO2+F_CH4*Cp_CH4);


%% Change vector converted to column vector
dzdM=dzdM';
end 
