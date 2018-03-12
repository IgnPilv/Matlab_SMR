% ----------------------------------------------------------------------- %
% calcz.m
% Cumene production packed bed reactor feed compressibility script
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
function z = calcz(Tz,Pz,Fb0,Fp0,Fc0,Fd0,Fi0,Ft0)
Tcrit=[562.16,364.76,631.15,689,369.82]; % K
Pcrit=[48.98,46.13,32.09,24.5,42.49]; % bar
Acent = [0.211,0.142,0.338,0.39,0.152]; % -
Ft0 = Fb0 + Fp0 + Fc0 + Fd0 + Fi0; % mol/s total inlet flow
[Tr]=Tz./Tcrit; % -
[Pr]=Pz./100000./Pcrit; % -
Trf = (Fb0/Ft0)*Tr(1) + (Fp0/Ft0)*Tr(2) + (Fc0/Ft0)*Tr(3) + (Fd0/Ft0)*Tr(4) + (Fi0/Ft0)*Tr(5);
% -
Prf = (Fb0/Ft0)*Pr(1) + (Fp0/Ft0)*Pr(2) + (Fc0/Ft0)*Pr(3) + (Fd0/Ft0)*Pr(4) + (Fi0/Ft0)*Pr(5);
% -
Acentf = (Fb0/Ft0)*Acent(1) + (Fp0/Ft0)*Acent(2) + (Fc0/Ft0)*Acent(3) + (Fd0/Ft0)*Acent(4) +(Fi0/Ft0)*Acent(5); % -
B1 = 0.083 - 0.422/(Trf^1.6); % -
B2 = 0.139 - 0.172/(Trf^4.2); % -
z = 1 + B1*(Prf/Trf) + Acentf*B2*(Prf/Trf); % Compressibility factor
end