% ----------------------------------------------------------------------- %
% MultipleReactorOptimizer.m
% Cumene production packed bed reactor design optimization script
% ----------------------------------------------------------------------- %
%
% Malcolm Bell, Sophie Fitchett, Adam Machon, Jack Palmer, Ignas Pilvinis
% S1433709 S1405595 S1426444 S1439172 S1452897
% 07/11/2017
% ----------------------------------------------------------------------- %
% SingleReactorPlotter.m
% MultipleReactorOptimizer.m
% calcz.m
% ODEsolver.m
% ----------------------------------------------------------------------- %
close all
clear all
clc
% -- DECLARATIONS OF INPUT PARAMETERS
global P0 T0 Ft0 rhoc U a z Cpcool mcool Rgas Rsi alpha
Wspan = linspace(0,20000,2000); % Range for the weight of the catalyst kg (xaxisof produced graphs)
Dt = 0.0789; % m tube diameter
% Nt = 700; % number of tubes (used for pressure drop)
P0 = 2500000; % Pa inlet pressure
T0 = 350+273; % K inlet temperature
Ta0 = 290+273; % K inlet coolant temperature
Fb0 = 70.8; % mol/s inlet benzene flow rate
Fp0 = 30.3; % mol/s inlet propylene flow rate
Fc0 = 0; % mol/s inlet cumene flow rate
Fd0 = 0; % mol/s inlet DIPB flow rate
Fi0 = 1.59474; % mol/s inlet propane flow rate
Ft0 = Fb0 + Fp0 + Fc0 + Fd0 + Fi0; % mol/s inlet total flow rate
rhoc = 1600; % kg/m³ particle density
% e = 0.5; % - voidage
U = 65; % W/(m² K) overall heat transfer coefficient
% Dp = 0.003; % m particle diameter
a = 4/Dt; % 1/m area to volume ratio
% Ac = 3.1415*Dt*Dt/4*Nt; % m² tube x-section area (used for pressuredrop)
% Cpcool = 2600; % J/(kg K)cooling heat capacity
% mcool = 25; % kg/s cooling mass flow
Rgas = 8.3144598; % J/molK
Rsi = 8.3144598; % kg/(m2 s2 K mol) gas constant in SI unitsfor density calculation
% -- PRESSURE DROP CALCULATIONS
% m0 = (Fb0*78.11+Fp0*42.08+Fc0*120.19+Fd0*162.276+Fi0*40.08)/1000; % kg/s Feed massflow (used for pressure drop)
% nu = 0.0000176939534465372; % kg/(m s)viscosity) Average value retrieved FROM THERMOWORKBENCH (used for pressure drop)
% rho0 = P0*(m0/Ft0)/(calcz(T0,P0,Fb0,Fp0,Fc0,Fd0,Fi0)*T0*Rsi); % kg/m³ (used forpressure drop)
% G = m0/(Ac); % kg/(m² s) massflux (through 1 tube) (used for pressure drop)
% Beta0 = (G/(rho0*Dp)) * ((1-e)/(e*e*e)) * (((150*(1-e)*nu)/Dp)+(1.75*G)); % (used forpressure drop)
% alpha = (2*Beta0)/(Ac*(1-e)*rhoc*P0); % (used forpressure drop)
% -- COMPRESSIBILITY CALCULATIONS
% z = calcz(T0,P0,Fb0,Fp0,Fc0,Fd0,Fi0);
% -- MULTIPLE INPUTS GENERATOR / SOLVER
k = 1; % dummy value for increments
T0 = 320+273; % C inlet temperature
for j = 1:1
Ta0 = 290+273; % C inlet coolant temperature
 for m0 = 1:1
 Fb0 = 30.3; % mol/s inlet benzene flowrate
 for l = 1:5
 z = calcz(T0,P0,Fb0,Fp0,Fc0,Fd0,Fi0); % - Compressibility factor
 i(:,k) = [Fb0; Fp0; Fc0; Fd0; Fi0; P0; T0; Ta0; 0;]; % Initial values for thedependent variables.
 [w(:,k),y(:,:,k)] = ode45(@ODEsolver,Wspan,i(:,k));
 Ti(k) = i(7,k)-273; % C inlet temperature
 Tai(k) = i(8,k)-273; % C inlet coolant temperature
 Pi(k) = i(6,k)/100000; % bar inlet pressure
 Fbi(k) = i(1,k); % mol/s inlet benzene flowrate
 Fpi(k) = i(2,k); % mol/s inlet propylene flowrate
 Wf(k) = interp1q(y(:,3,k),w(:,k),28.93); % kg mass of catalyst
 Tf(k) = interp1(w(:,k),y(:,7,k),Wf(k))-273; % C output temperature
 Taf(k) = interp1(w(:,k),y(:,8,k),Wf(k))-273; % C output coolant temperature
 Pf(k) = interp1(w(:,k),y(:,6,k),Wf(k))/100000; % bar output pressure
 Lf(k) = Wf(k)/(Ac*(1-e)*rhoc); % m length of the reeactor
 conv(k) = interp1(w(:,k),y(:,3,k),Wf(k))/i(2,k); % - conversion of propylene tocumene
 Fdf(k) = interp1(w(:,k),y(:,4,k),Wf(k)); % mol/s output DIPB flow
 Vf(k) = Wf(k)/(rhoc*(1-e)); % m³ reactor volume
 dPbar(k) = Pi(k)-Pf(k); % bar pressure drop
 Fb0 = Fb0 + 10;
 k=k+1
 end
 Ta0 = Ta0 + 10;
 end
T0 = T0 + 10;
end
% -- TABULAR DATA
T = table(Ti.',Tai.',Pi.',Fbi.',Fpi.',Wf.',Tf.',Taf.',Pf.',Lf.',conv.',Fdf.',Vf.',dPbar.');
T.Properties.VariableNames = {'T0' 'Ta0' 'P0' 'Fb0' 'Fp0' 'W' 'Tf' 'Taf' 'Pf' 'L' 'X' 'Fdf' 'Vf' 'dPf'}
% -- DATA EXPORTER
writetable(T,'myData.xls');