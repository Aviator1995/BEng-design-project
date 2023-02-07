% Preliinary thrust calculation, Straight & Level, Unaccelerated Flight.
% For altitude at Troposphere & Tropopause (i.e. Max ceiling 25Km)
clear all; clear; clc

%%%%%%%%%%%%%%%%% Initial parameters
Cd     = 0.014;  % Drag coefficient at flight condition
m      = 50000;  % Mass of aircraft (Kg)
W      = m*9.81; % Weight of aircraft(N) (T/W normalised under this value)
M      = 1.8;    % Mach number
S      = 82.9;   % Wing area (m2)
H      = 0;   % Altitude (km)
T_pe_W = 0.17;   % Thrust (per engine) / Weight   <== Market research

%%%%%%%%%%%%%%%%%% Thermodynamic properties
gamma = 1.4;     % Speific heat ratio
Rg    = 287;     % gas constant(J/Kg.K)
rho_0 = 1.225;   % Density at sea level (Kg/m3)
T_sls = 15.04;   % Temparature at sea level (degree Cel.)
P_sls = 101.325; % Pressure sea level (kPa)

%% Section 1 - Altitude conditions

h      = H*1000; % Altitude (m)

if H <= 11, % Atmosphere - Troposhpere 
    
    T_alt = 288.15-0.0065.*h;      % Temp. (K) in Troposphere
    fprintf('Temparature at altitude (K): %.3f \n',T_alt)
    theta = 1-0.000022558.*h; % density lapse ratio in troposphere 
    del = (1-0.000022558.*h).^5.2558;% Hydrostatic equation (pressure var.)
    P_alt = P_sls*del;   % Pressure at altitude (kPa)
    fprintf('Pressure at altitude (kPa): %.3f \n',P_alt)
    sig = del/theta; % density ratio
    rho = sig*rho_0; % density at altitude (kg/m3)
    fprintf('Density at altitude (Kg/m3): %.3f \n',rho)
    
elseif H > 11, % Atmosphere - Tropopause 
    T_alt = 216.65;            % Temp. (K) in tropopause 
    fprintf('Temparature at altitude (K): %.3f \n',T_alt)
    theta = 0.7519;         % density lapse ratio in tropopause 
    del = .22336.*exp(1.7346-0.00015769.*h);% pressure variation
    P_alt = P_sls*del;   % Pressure at altitude (kPa)
    fprintf('Pressure at altitude (kPa): %.3f \n',P_alt)
    sig = del/theta; % density ratio
    rho = sig*rho_0; % density at altitude (kg/m3)
    fprintf('Density at altitude (Kg/m3): %.3f \n',rho)
end

% Flight Velocity
c     = sqrt(gamma.*Rg.*T_alt);  % speed of sound (m/s)
V     = M.*c;                    % Flight velocity (m/s)
fprintf('Aircraft Velocity (m/s): %.3f \n',V)

% Section 2 - Thrust Required
D     = (1/2).*rho.*V^2.*S.*Cd;  % Drag equation
T_req = D;                       % S&L unaccelerated
fprintf('Drag/Static Thrust required at flight condition (N): %.3f \n',T_req)
% **Note, Flight condition at max vlocity S&L

% Section 3 - T/W
T_W = T_req./W;                  % Thrust to Weight 
fprintf('Thrust to Weight Ratio: %.3f \n',T_W)

% Section 4 - Thrust per engine
T_pe = T_pe_W.*W;                % Thrust per engine (N)
fprintf('Thrust per engine (N): %.3f \n',T_pe)

% Section 5 - Number of engines
N = T_req/T_pe;                  % Number of engines required
fprintf('Number of engines: %.3f \n',N)
N_act = ceil(N); % Actual number of engines
fprintf('Actual Number of engines: %.3f \n',N_act)

fprintf('Altitude (km): %.3f \n',H)
