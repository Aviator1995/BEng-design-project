%% Take off Parammeters and Distance
clc
clear all

%% Input parameters

g     = 9.81;         % gravity (m/s2)
Wto   = 44000*g;      % Weight of drone (N)
S     = 90;           % wing area (m2)
b     = 6.46.*2;      % span of wings (m)
AR    = 1.86;         % Aspect Ratio
rho   = 1.225;        % density (kg/m3)
CLmax = 0.9;          % coefficient of lift
Cd0   = 0.005;        % coefficient of zero-lift drag
aoa   = 3;            % Angle of Attack
h     = 5;            % height of wings from the ground
c     = 13.94;        % Length of chord
mu    = 0.05;         % friction for concrete
T     = 222000;       % Max Thrust (N)
Ho    = 10;           % height of obstacle
CL_ge = 0.74;         % Coefficient of Lift at take-off
Vin   = 0;            % Initial velocity
e     = 0.7;          % oswalds coefficient
t_r   = 3;            % rotation time(s)

%% Take-off Velocity

% Take-off velocity
Vto =1.2*(sqrt((2*Wto)/(S*rho*CLmax)));

% Lift at Ground Roll w/ ground effects
Lg= (CL_ge*rho*(Vto.^2)*S*0.5);      % Lift at Ground Roll

% Drag at Ground Roll (considering ground effect)
phi=((16*(h/b))^2)/(1+(16*(h/b))^2); % Ground Effect Coefficient
Cd=Cd0+((phi*(CL_ge^2))/(pi*e*AR));  % Drag Coefficient
D=0.5*rho*S*Cd*Vto.^2;               % Drag, N

% Phase 1 - Ground Roll

F   = mu*(Wto-Lg);              % Friction force, N
a   =(g/Wto)*(T-D-F);           % acceleration (m/s2)
SG  = ((Vto^2)/(2*a)) - ((Vin^2)/(2*a));    % Ground roll
fprintf('Ground Roll Distance (m): %.3f \n',SG)

% Phase 2 - Rotation

SR=t_r*Vto;           % Horizontal distance during rotation (SR)
fprintf('Horizontal rotation distance (m): %.3f \n',SR)

% Phase 3 - Transition

gamma= asind((T-D)/Wto);
gammaCl = gamma;            % Climb angle
Rtr=((Vto^2)/(0.15*g));
Str=Rtr*sind(gammaCl);
fprintf('Transition distance (m): %.3f \n',Str)

% Phase 4 - Climb

Htr=Rtr*(1-cosd(gammaCl));
SCl= ((Ho-Htr)/(tand(gammaCl)));
fprintf('Climb distance (m): %.3f \n',SCl)

% Phase 5 - Total take-off distance

Sto=SG+SR+Str+SCl;
fprintf('Total take off distance (m): %.3f \n',Sto)

%% Plot Results

% Phase 1 plot
SG_x = [0,SG]; 
SG_y = [0,0];
plot(SG_x, SG_y,'DisplayName','Ground Roll','LineWidth',1.5); hold on;     

% Phase 2 plot
SR_x = [SG, (SR+SG)]; 
SR_y = [0,0];
plot(SR_x, SR_y,'DisplayName','Rotation','LineWidth',1.5); hold on;     

% Phase 3 plot
r = Rtr; b = r;
n = 50;                   % No. distance interval transition stage
y = zeros(n+1,1);         % height vector

n_ini = SR+SG;                        % initial distance, transition
n_inc = (((Str+SR+SG)-(SR+SG))./n);   % distance increment
n_final = (SG+SR+Str);                % fina distance, transition

x = n_ini:n_inc:n_final; x = x'; % horizontal distance
a = SR+SG;

% corresponding y value
for i = 1:n+1, y(i) = sqrt((r.^2)-((x(i) - a).^2)) + b;end

y = y-b; 
y = r-y;     % Invert the y vector (opposite circ)
Str_x = x; 
Str_y = y;    
plot(Str_x, Str_y, 'k', 'DisplayName','Transistion','LineWidth',1.5); hold on;

% Phase 4 plot
SCL_x = [Str+SR+SG, Sto]; SCL_y = [y(n+1), Ho];
plot(SCL_x, SCL_y,'DisplayName','Climb','LineWidth',1.5);legend('show');
xlabel('Distance (m)'); ylabel('Altitude (m)');
grid on; axis equal; hold off
title('Take-Off flight path')
