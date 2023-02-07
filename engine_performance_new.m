clear all; clear; clc;

TA_0     = 222000;         % Thrust available at sea-level (N)
S        = 90;             % Wing area (m2)
l        = 35;             % length
CD0      = 0.005;          % Drag Coefficient, zero lift
M_max    = 3;              % Max. Mach No.
v_min    = 50;             % Stall velocity (m/s)
n        = 1;              % load factor
m        = 44000;          % mass of Aircraft (kg)
g        = 9.81;           % gravitational constant (m/s2)
W        = m * g;          % Weight of Aircraft (N)
AR       = 1.86;           % Aspect Ratio
H_max    = 15.24;          % Cruise Ceiling (Km) <-- 50,000 ft
H_sls    = 0;              % Gound level
rho_0    = 1.225;          % Air Density (Kg/m3) 
AR       = 1.86;           % Aspect Ratio
lambda   = 65;             % .25c Sweep angle 
b        = 6.46.*2;        % Wing span (m)
Vol      = 36;             % Wing volume (m3)

%% Section 1 - Power available and power required (Max velocity v altitude)

% Power available at sea level v velocity
M = M_max;
[ rho, T_alt ] = altitude_properties( H_sls );
[ V_max,c ]      = Mach_Velocity( T_alt, M );

fprintf('Temparature at Sea Level (K): %.3f \n',T_alt);
fprintf('Density at Sea Level (Kg/m3): %.3f \n',rho);

i = 1; V_max = ceil(V_max);
for v = 0:V_max,
    V_sls_0(i) = v;
    Mach_sls_0(i) = V_sls_0(i)./c;
    PA_0(i) = TA_0 .* v;
    i = i+1;
end

% Power required at sea level
[ Mach_sls, D_sls, V_sls ] = Supersonic_drag_approx( H_sls, l, S, W, CD0, v_min, M, b, Vol);
Pr_0 = D_sls .* V_sls;

% Create Thrust available at sea level vector
j = length(Mach_sls); for jj=1:j, T_A_0(jj)=TA_0;end

%% Display Results - Power available/required at Sea level

figure(1), 

subplot(1,2,1);  % Plot Power available/required v Mach at sea level
plot( Mach_sls_0, PA_0, 'k' ); hold on   % Plot Power available v Mach
subplot(1,2,1); grid minor
plot( Mach_sls, Pr_0 , '--k'); hold off  % Plot power required  v Mach
title('Power v Mach, at Sea level'); ylabel('Power (W)');xlabel('Mach')
legend( 'P_A_o' , 'P_R_o' , 'location' , 'Northeast' );

subplot(1,2,2);  % Plot Power available/required v Velocity
plot( V_sls_0, PA_0, 'k' );  hold on    % Plot Power available v Velocity
subplot(1,2,2); grid minor
plot( V_sls, Pr_0, '--k' );  hold off   % Plot Power required  v Velocity
title('Power v Velocity, at Sea level '); ylabel('Power (W)');
xlabel('Velocity (m/s)'); legend('P_A_o','P_R_o','location','Northeast');
hold off

figure(2),

subplot(1,2,1);  % Plot Thrust available/required v Mach Number 
plot(Mach_sls, T_A_0, 'k'); hold on     % Plot Thrust available v Mach 
subplot(1,2,1); grid minor
plot(Mach_sls, D_sls, '--k'); hold off  % Plot Thrust required  v Mach
title('Thrust v Mach, at Sea level'); 
xlabel('Mach'); ylabel('Thrust (N)');
legend('T_Ao' , 'T_Ro', 'location' , 'Northeast');

subplot(1,2,2);  % Plot Thrust available/required v Velocity
plot(V_sls, T_A_0, 'k'); hold on        % Plot Thrust available v Velocity
subplot(1,2,2); grid minor
plot(V_sls, D_sls, '--k'); hold off     % Plot Thrust required  v Velocity 
title('Thrust v Velocity'); ylabel('Thrust (N)'); xlabel('Velocty (m/s)');
legend('T_Ao' , 'T_Ro', 'location' , 'Northeast');
hold off

%% Power available at altitude v velocity
[ rho, T_alt ] = altitude_properties( H_max );
[ V_max,c ]      = Mach_Velocity( T_alt, M );

fprintf('Temparature at altitude (K): %.3f \n',T_alt);
fprintf('Density at altitude (Kg/m3): %.3f \n',rho);

TA_alt = (rho./rho_0) .* TA_0; % Thrust available at altitude

i = 1; V_max = ceil(V_max);
for v = 0:V_max,
    V_alt(i) = v;
    Mach_alt_r(i) = V_alt(i)./c;
    PA_alt(i) = TA_alt .* v;
    i = i+1;
end

% Power required at altitude
[ Mach_alt, D_alt, V_alt_r ] = Supersonic_drag_approx( H_max, l, S, W, CD0, v_min, M, b, Vol);
Pr_alt = D_alt .* V_alt_r;

% Create Thrust available at altitude vector
j = length(Mach_alt); for jj=1:j, T_A_alt(jj)=TA_alt;end

%% Display Results - Power available/required at Altitude
format shortG
H_max_ft = H_max.*3280.84;

figure(3), 

subplot(1,2,1); % Plot Power available/required v Mach at altitude
plot( Mach_alt_r, PA_alt, 'k' ); hold on % Plot power available v Mach
subplot(1,2,1); grid minor
plot(Mach_alt, Pr_alt, '--k'); hold off  % Plot power required  v Mach
title(['Power v Mach, at ' num2str(H_max_ft) ' feet']); 
xlabel('Mach'); ylabel('Power (W)');
legend( 'P_A' , 'P_R' , 'location' , 'Northeast' );

subplot(1,2,2); % Plot Power available/required v Velocity at altitude
plot(V_alt, PA_alt, 'k'); hold on        % PLot power available v Velocity
subplot(1,2,2); grid minor
plot(V_alt_r, Pr_alt, '--k' ); hold off  % PLot power required  v Velocity 
title(['Power v Velocity, at ' num2str(H_max_ft) ' feet']);
xlabel('Velocity (m/s)'); ylabel('Power (W)');
legend( 'P_A' , 'P_R' , 'location' , 'Northeast' );

figure(4),

subplot(1,2,1); % Plot Thrust available/required v Mach at altitude
plot(Mach_alt,T_A_alt, 'k'); hold on    % Plot Thrust available v Mach 
subplot(1,2,1); grid minor
plot(Mach_alt, D_alt, '--k'); hold off  % Plot Thrust required  v Mach 
title(['Thrust v Mach, at ' num2str(H_max_ft) ' feet']);
ylabel(' Thrust (N)'); xlabel(' Mach '); ylim([0 TA_0]);
legend('T_A', 'T_R' , 'location', 'Northeast');

subplot(1,2,2); % Plot Thrust available/required v Velocity at altitude
plot(V_alt_r,T_A_alt, 'k'); hold on     % Plot Thrust available v Velocity 
subplot(1,2,2); grid minor
plot(V_alt_r,D_alt, '--k'); hold off    % Plot Thrust required  v Velocity 
title(['Thrust v Velocity, at ' num2str(H_max_ft) ' feet']);
ylabel('Thrust (N)'); xlabel('Velocity (m/s)');ylim([0 TA_0]);
legend('T_A', 'T_R' , 'location', 'Northeast');
hold off

%% ROC - v : ROC == (PA - PR)/ MTOW

i=1;
h_max = H_max+5;
for h = H_sls:0.01:h_max;
    
    [ Mach, D, V ] = Supersonic_drag_approx( h, l, S, W, CD0, v_min, M, b, Vol);
    
    Power_r = D.*V;
    
    [ rho, T_alt ] = altitude_properties( h );
    
    Thrust_A = (rho./rho_0) .* TA_0; % Thrust available at altitude
    Power_A = Thrust_A.*V;
    
    del_P = Power_A - Power_r;
    del_P_max(i) = max(del_P);
        
    i=i+1;
end
ROC_max = del_P_max./(W*10);

figure(5),
plot(ROC_max, H_sls:0.01:h_max, 'k'); title('Maximum ROC v Altitude');
xlabel('Rate of Climb (m/s)'); ylabel('Altitude (km)');
xlim([0 max(ROC_max)]);
grid minor; hold off


