clc
clear all;

% MatLab script for modelling of a PSA system for purification of Hydrogen
% containing CO2, CO, and CH4 impurities. 
%
% Assumptions
%
% The CO2 and CH4 are both adsorbed on AC and constitute 90% of the
% impurities. Assume that the whole impurity acts the same way. 
% To account for the 10% of CO, a layer of zeolite 5A can be added at the
% top of the column

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants used in the calculation
K_A = 1.578;
K_B = 0.0386;
epsilon = 0.3;
rho_b = 500; % kg./m3
y_F = 0.056;
mu = 4.92e-5; % Pas
M_F = 0.0022; % kg./mol

% Since the isotherms are linear, beta = theta
beta_A = 1./(1 + ((1 - epsilon).* K_A)./epsilon);
beta_B = 1./(1 + ((1 - epsilon).* K_B)./epsilon);
beta = beta_A./beta_B;

theta_A = beta_A;
theta_B = beta_B;
theta = theta_A./theta_B;

% Molar feed flowrate and density 
rho_F = 0.89; % kmol./m3 (see excel 'Cycle Steps')
Q = 540; % mol./s

% Pressures, temperature and gas constant
P_H = 21; % atm
P_L = 1; % atm
P_LPR = 17.25; % atm
R = 8.2060e-5; % m3atm./molK
T = 293.15; % K

% phi parameter and Reynolds number
APR = P_H./P_L; 
phi_tF = Q./APR; % mol./s
phi_Vads = (epsilon.*P_L)./(beta_A.*R.*T); % mol./m3
Vads_tF = phi_tF./phi_Vads; % m3./s

d_p = 0.0025; % m 
Re = [12.9; 20; 30; 50; 1000];
A_c = (Q.*M_F.*d_p)./(mu.*Re); % m2

L_tF = Vads_tF./A_c; % m./s
t_F = 100; %s
L = L_tF.*t_F; % m
d = sqrt((4.*A_c)./pi); % m
d_industrial = d./0.3048; % m
CYL_RATIO = L./d; 

% Recovery
Q_in_tF = (epsilon.*A_c.*L.*P_H)./(theta_A.*R.*T); % mol
Q_in_tPR = ((epsilon.*A_c.*L)./(beta_A.*R.*T)).*(P_H-P_LPR); % mol
Q_out_tF = (1+(theta-1).*y_F).*epsilon.*A_c.*(P_H./(R.*T)).*(L./theta_A); % mol

% no purge in recovery calc - the product is not split for purge only for
% pressurization stage of another column

% Recovery
R_B = (Q_out_tF - Q_in_tPR)./(Q_out_tF.*(1-y_F));

% Find the mass and volume of adsorbent needed to purify the hydrogen (10
% beds needed)
V = 0.25.*pi.*d.^2.*L; % m3
mass = rho_b.*V; % kg

mass_total = 10.*mass; % kg

% Adjusting the dimensions of the column for inert packing and the zeolite
% 5A layer

epsilon_5A = 0.35; 
rho_5A = 1180; % kg./m3
ratio_AC5A = 0.25./0.75;
mass_5A = mass.*(ratio_AC5A); % kg
total_mass_5A = 10.*mass_5A; % kg
volume_5A = mass_5A./rho_5A; % m3
volume_5A_wvoid = volume_5A./(1-epsilon_5A); % m3
V_inert = (V+volume_5A_wvoid)./0.9; % m3

Volume_total = V_inert; % m3
L_adjusted = Volume_total./(A_c); % m
Adj_cyl_ratio = L_adjusted./d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical design
% Design for both a torisphere and an ellipsoidal head to see which one is
% more economical

design_P = (21 - 1).*1.1; % atm (take as 10% above operating gauge pressure)
% Convert atm to bar
design_P_bar = design_P .* 1.01325; % bar
% Convert to N./mm2
design_P_N = design_P_bar .* 0.1; % N./mm2
% Using carbon steel, the maximum allowed stress at 100 F (37.8 C) 
stress = 88.9; % N./mm2

% Cylindrical section
Eff = 1; % weld efficiency
thickness_c = (design_P_N.*d.*10.^3)./(2.*stress.*Eff-1.2.*design_P_N)+2; % mm (+2 mm to allow for corrosion). This could be rounded up to 1 inch

% Domed section
% Dished head (torosphere)
R_C = d; % m
knuckle_radius = 0.06.*R_C; % m
thickness_t = (0.885.*design_P_N.*R_C.*10..^3)./(stress.*Eff - 0.1.*design_P_N); % mm

% Ellpisoidal head (ratio of major:mior axes is 2:1)
thickness_e = (design_P_N.*d.*10..^3)./(2.*stress.*Eff - 0.2.*design_P_N); % mm

% The column should therefore use an ellipsoidal head with a uniform
% thickness of plate of 1 inch across both sections

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displayed values

fprintf('The column diameter is %s m \n', d)
fprintf('The total column height is %s m \n', L_adjusted)
fprintf('The mass of AC needed is %s kg \n', mass_total)
fprintf('The mass of zeolite 5A required is %s kg \n', total_mass_5A)
fprintf('Volume of inert packing is %s m3 \n', V_inert - volume_5A_wvoid - V)
fprintf('Reynolds number is %i because of the implication \n', Re)
fprintf('Hydrogen recovery is %s percent \n', R_B.*100)
fprintf('Cylinder ratio is %s \n', Adj_cyl_ratio)
disp(d);
disp(Re);
plot(Re,L_adjusted);


