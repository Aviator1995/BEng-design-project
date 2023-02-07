function [ rho, T_alt ] = altitude_properties( H )
% Computing atmospheric air density and temparature for sls & altitude in
% the Troposphere & Tropopause (i.e. Max ceiling 25Km)

%% ALTITUDE PROPERTIES
%%%%%%%%%%%%%%%%%% Thermodynamic properties
gamma = 1.4;     % Speific heat ratio
Rg    = 287;     % gas constant(J/Kg.K)
rho_0 = 1.225;   % Density at sea level (Kg/m3)
T_sls = 15.04;   % Temparature at sea level (degree Cel.)
P_sls = 101.325; % Pressure sea level (kPa)

%Altitude conditions

h      = H*1000; % Altitude (m)

if H <= 11, % Atmosphere - Troposhpere 
    
    T_alt = 288.15-0.0065.*h;      % Temp. (K) in Troposphere
    theta = 1-0.000022558.*h; % density lapse ratio in troposphere 
    del = (1-0.000022558.*h).^5.2558;% Hydrostatic equation (pressure var.)
    P_alt = P_sls*del;   % Pressure at altitude (kPa)
    sig = del/theta; % density ratio
    rho = sig*rho_0; % density at altitude (kg/m3)
    
elseif H > 11, % Atmosphere - Tropopause 
    T_alt = 216.65;            % Temp. (K) in tropopause 
    theta = 0.7519;         % density lapse ratio in tropopause 
    del = .22336.*exp(1.7346-0.00015769.*h);% pressure variation
    P_alt = P_sls*del;   % Pressure at altitude (kPa)
    sig = del/theta; % density ratio
    rho = sig*rho_0; % density at altitude (kg/m3)
end

end

