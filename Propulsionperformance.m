% Code for computing T-v, p-V, ROC-V, ROC-h.
clear all; clear; clc

disp(['*****************************************************************'])
disp(['*                Propulsion, parameter sweep                    *'])
disp(['* ************************************************************* *'])
disp(['*   Section 1   - Initial constants & parameters                *'])
disp(['*   Section 1.1 - Prompt, parameter questions & CL,CD print     *'])
disp(['*   Section 2   - P-v diagram                                   *'])
disp(['*   Section 3   - T-v diagram                                   *']) 
disp(['*   Section 4   - ROC - v diagram                               *'])
disp(['*   Section 5   - ROC - Altitude, h  diagram                    *'])
disp(['* ************************************************************* *'])
disp(['* ------------------------------------------------------------- *'])
disp(['*   *All values in SI units                                     *'])
disp(['*   Enter parameters from section 1 [0] OR section 1.1 [1]      *'])
disp(['*                                                               *'])
disp(['*****************************************************************'])
disp(['                                                                 '])
disp('===================================================================')
disp('Parameter input')
disp('                                                                   ')
disp('From Script, Manual input: [0]')
disp('Parameter input from Menu: [1]')
disp('                                                                   ')

% Enter option of parameter input method
prompt = 'Method of parameter input [0|1]: '; option = input(prompt);

% 0 - Section 1   (Script)
% 1 - Section 1.1 (Questions)

disp(['                                                                 '])

%% Section 1 - Initial constants & parameters

if option == 0;

    PA       = 163323624.421;  % Power available (W)
    PA_0     = 190773017.588 ; % Power available at sea-level (W)
    rho_0    = 1.225;          % Air density  at sea-level(kg/m3)
    S        = 82.9;           % Wing area (m2)
    CD_0     = 0.01;           % Drag Coefficient, zero lift
    n        = 1;              % load factor
    m        = 50000;          % mass of aircraft (kg)
    g        = 9.81;           % gravitational constant (m/s2)
    e        = 0.7;            % Oswald eff. factor
    AR       = 7;              % Aspect Ratio
    V_min    = 80;             % Min. velocity (m/s)
    V_inc    = 1;              % Veocity increment (m/s)
    V_max    = 617.4;          % Max. Velocity (m/s)
    V_cruise = 583.1;          % Velocity_Cruise (m/s)
    h_min    = 0;              % Minimum altitude (Km)
    h_max    = 18.3;           % Maximum altitude (Km)

% Section 1.1 Prompt

elseif option == 1;
    
disp(['"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""'])
disp(['             Enter values for following parameters               '])
disp(['-----------------------------------------------------------------'])
disp(['                                                                 '])

    prompt = 'Power available, PA (W): ';               PA = input(prompt);
    prompt = 'Power available at sea-level, PA_0 (W): ';PA_0=input(prompt);
    prompt = 'Air density at sea-level (kg/m3): ';   rho_0 = input(prompt);
    prompt = 'Wing area (m2): ';                         S = input(prompt);
    prompt = 'Drag Coefficient - zero lift, CD_0: ';  CD_0 = input(prompt);
    prompt = 'load factor, n: ';                         n = input(prompt); 
    prompt = 'mass of aircraft, m (kg): ';               m = input(prompt);
    prompt = 'gravitational constant (m/s2): ';          g = input(prompt);
    prompt = 'Oswald efficiency factor, e: ';            e = input(prompt);
    prompt = 'Aspect ratio: ';                          AR = input(prompt);
    prompt = 'Minimum velocity (m/s): ';             V_min = input(prompt);
    prompt = 'Velocity increment (m/s): ';           V_inc = input(prompt);
    prompt = 'Maximum velocity (m/s): ';             V_max = input(prompt);
    prompt = 'Cruise speed (m/s): ';              V_cruise = input(prompt);
    prompt = 'Minimum Altitude (km): ';              h_min = input(prompt);
    prompt = 'Maximum Altitude (km): ';              h_max = input(prompt);

end

% Resultant parameters

format longG

CL = sqrt(3*pi*e*AR*CD_0);         % Lift Coefficient
CD = CD_0+((CL.^2)/(pi.*e.*AR));   % Drag Coefficient
W = m.*g;                          % Weight of Aircraft (N)

% Print resultant parameters

disp(['                                                                 '])
disp(['"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""'])
disp(['-----------------------------------------------------------------'])
disp(['Output parameters'])
disp(['                                                                 '])
1

fprintf('Lift Coefficient, CL: %.3f \n',CL);
fprintf('Drag Coefficient, CD: %.3f \n',CD);
fprintf('Aircraft Weight (N): %.3f \n',W);

%% Section 2 P-v diagram

V = V_min:V_inc:V_max; % Velocity range (m/s)

% Thrust required for Vel.
Tr = (0.5.*rho_0.*S.*CD_0.*(V.^2))+(((2.*((n.*W).^2))./...
    (pi.*e.*AR.*rho_0.*S)).* (1./(V.^2)));

Tr_min = min(Tr);

% Power required

Pr = Tr .* V;

Pr_min = min(Pr);

% At cruise

Tr_c = (0.5.*rho_0.*S.*CD_0.*(V_cruise.^2))+(((2.*((n.*W).^2))./...
    (pi.*e.*AR.*rho_0.*S)).* (1./(V_cruise.^2)));     % Tr at cruise flight

Pr_c = Tr_c .*V_cruise;   % Pr at cruise flight

% Plot P-V

figure(1); hold on;
plot(V, Pr, 'k'); % plot full range graph
hold on

xlabel('Velocity (m/s)');
ylabel('Power Required, Pr (W)');
title('P-v Diagram');
grid on
hold on

plot(V_cruise, Pr_c, '*r'); hold on;

plot([V_min, V_cruise], [Pr_c, Pr_c], '--r'); hold on;

Pr_ylabel = Pr_min-20; Pr_ylabel = round(Pr_ylabel, -1);
plot([V_cruise, V_cruise], [Pr_ylabel, Pr_c], '--r'); hold on;

ylim([Pr_ylabel, inf]); hold off;

% Print

Pr_c = round(Pr_c, 3);   % Print power required at cruise

disp(['                                                                 '])
disp(['        -------------------------------------------------        '])
disp(['                                                                 '])

fprintf('SECTION 2 - Power required at cruise flight (W): %.3f \n',Pr_c)

%% Section 3 T-v diagram 

% Plot T-V

figure(2); hold on;
plot(V, Tr, 'k'); hold on; % Plot full range graph

xlabel('Velocity (m/s)');
ylabel('Thrust Required, Tr (N)');
title('T-v Diagram');
grid on
hold on

plot(V_cruise, Tr_c, '*r'); hold on;

plot([V_min, V_cruise], [Tr_c, Tr_c], '--r'); hold on;

Tr_ylabel = Tr_min-2; Tr_ylabel = round(Tr_ylabel, -1);
plot([V_cruise, V_cruise], [Tr_ylabel, Tr_c], '--r'); hold on;

ylim([Tr_ylabel, inf]); hold off;

% Print

Tr_c = round(Tr_c, 3);  % Print thrust required at cruise
fprintf('SECTION 3 - Thrust required at cruise flight (N): %.3f \n',Tr_c)



%% ------------------------------------------------------------------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%   ROC - V,h diagram   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% -------------------------------------------------------------------------

%% Section 4 - ROC - V diagram

ROC = (PA - Pr)./W;

figure(3); hold on;
plot(V, ROC, 'k'); grid on; hold on;
xlabel('Velocity (m/s)'); ylabel('ROC (m/s)');
title('ROC - v'); hold on;

% ROC max

ROC_max = max(ROC);
P_T = [Pr; Tr];
P_T = P_T'; % MATRIX ARRANGE AS P_T

Tr_p = P_T(P_T(:,:)==Pr_min,:);
Tr_p = Tr_p(:,2);

v_ROCmax = Pr_min/Tr_p;

plot(v_ROCmax, ROC_max, '*r'); hold off;

fprintf('SECTION 4 - ROC max (m/s): %.3f \n',ROC_max) 
fprintf('SECTION 4 - Velocity at ROC_max (m/s): %.3f \n',v_ROCmax)

%% Section 5

% ROC - Altitude, h  Diagram

h = h_min:0.1:h_max;
PA = ((20-h)./(20+h)).*PA_0;
rho = ((20-h)./(20+h)).*rho_0;
PR = sqrt((2.*(W.^3))./(rho.*S)).*(CD/CL.^1.5);
ROC_max = (PA-PR)./W;

figure(4); hold on;
plot(ROC_max,h,'k'); grid on;
xlabel('ROC max (m/s)'); ylabel('Altitude, h (km)');
title('ROCmax(h)');
legend('\it ROCmax(h)','location','Northeast');
xlim([0, inf]); hold off;


disp(['                                                                 '])
disp(['"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""'])
disp(['                                                                 '])
disp('Quit:           [0]')
disp('Return to Menu: [1]')
disp('                                                                   ')
prompt = 'Back to menu: '; option = input(prompt);

if option == 1; clear; run('Propulsionperformance.m'); end
if option == 0; clear; clc; end
