clc
clear all

n = 1;
W = 120;
rho = 1.225;
S = 0.624;
CL_max = 1.55;

V = 10:0.1:25;
% Equation of lift to determine the Value of Vs for given input parameters
%   Input arg: 
% CL - coefficient of lift
% n - Load factor
% W - Aircraft weight (kg)
% rho - air density (kg/m3)
% V - speed (m/s)
% S - wing surface area (m2)


figure(1),

i = 1;

for V_i = 10:0.1:25;
    CL(i) = (2.*n.*W)./(rho.*S.*(V_i.^2));
    i=i+1;
end

plot(V,CL); hold on

V_s = sqrt((2.*n.*W)./(rho.*S.*CL_max));
plot(V_s, CL_max, '*r'); hold on; grid minor
xlabel('Airspeed (m/s)'); ylabel('Lft cOEFF'); title('Landing speed (V_TD)');


