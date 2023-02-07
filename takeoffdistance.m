%% Take off Parammeters and distance
clear all; clear; clc

% Input parameters

a  = 20; %Acceleration
W  = 15; % Weight of drone
S  = 1; %wing area
rho=1.225; %density
CLmax=1.5; %coefficient of lift
g=9.81; %gravity
mu=.6; %friction
Wto=14.99; %take off weight
T=15; %thrust (N)
D=6; %drag
Lg=12; %lift at ground roll
Ho=80; %height of obstacle
V_ini = 0; % Initial velocity

%% Stage 1 - Ground Roll

Vto = 1.2*((2*W)/(S*rho*CLmax)); % Take-off velocity
% V=Vto; %Velocity
F   = mu*(Wto-Lg); % Friction force
a   =(g/Wto)*(T-D-F); % acceleration
SG  = ((Vto^2)/(2*a)) - ((V_ini^2)/(2*a)); % Ground roll
fprintf('Ground Roll (m): %.3f \n',SG)

%% Stage 2 - Rotation

SR=3*Vto; % Horizontal distance during rotation (SR)
fprintf('Horizontal rotation distance (m): %.3f \n',SR)

%% Stage 3 - Transition

gamma=asind((T-D)/W); gammaCl = gamma; % Climb angle
Rtr=((Vto^2)/(0.15*g));
Str=Rtr*sind(gammaCl);
fprintf('Transition distance (m): %.3f \n',Str)

%% Stage 4 - Climb
Htr=Rtr*(1-cosd(gammaCl));
SCl= ((Ho-Htr)/(tand(gammaCl)));
fprintf('Climb distance (m): %.3f \n',SCl)


%% Plot Results

Sto=SG+SR+Str+SCl;
fprintf('Total take off distance (m): %.3f \n',Sto)

% Stage 1 plot
SG_x = [0,SG]; SG_y = [0,0];
plot(SG_x, SG_y,'DisplayName','SG'); hold on;     % plot SG

 % Stage 2 plot
SR_x = [SG, (SR+SG)]; SR_y = [0,0];
plot(SR_x, SR_y,'DisplayName','SR'); hold on;     % plot SR

% Stage 3 plot
r = Rtr; b = r;
n = 50; % Number of node for transition stage
y = zeros(n+1,1);% height array
x = SR+SG:(((Str+SR+SG)-(SR+SG))./n):(SG+SR+Str); x = x'; % horizontal distance
a = SR+SG;

for i = 1:n+1,
    y(i) = sqrt((r.^2)-((x(i) - a).^2)) + b; % calc corresponding y value
end

y = y-b; y = r-y; % Invert the y values
Str_x = x; Str_y = y;    
plot(Str_x, Str_y, 'k', 'DisplayName','Str'); hold on; % Plot STR

% Stage 4 plot
SCL_x = [Str+SR+SG, Sto]; SCL_y = [y(n+1), Ho];
plot(SCL_x, SCL_y,'DisplayName','SCl');legend('show');
xlabel('Distance (m)'); ylabel('Altitude (m)');
grid on; axis equal; hold off
