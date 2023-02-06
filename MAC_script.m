clear all
close all
clc

AR = 13;
c = 0.225;
t = 0.3;
theta = 0.5; % Radians

%% Determiing planform coordinates;
t_c = c.*t;  % Tip cord of main wing
b = (AR.*(c.*(1+t)))./2; % Main wing span.
h = b./2; % Main wing half-span
plot(0,0,'pr');hold on;

sweep = h.*tan(theta); % sweep tip, quarter cord (y-axis)
c_4 = c.*0.75; % quarte cord from aft of aerofoil (y-axis)

ct_4 = c_4 - sweep; % tip quarter cord (y-axis)
ct_l = ct_4 + (t_c.*0.25); % tip cord leading edge
ct_t = ct_4 - (t_c.*0.75); % tip cord trailing edge

% Root cord
Root_y = [c, c_4, 0];
Root_x = zeros(1,3); 
%Tip cord
Tip_y = [ct_l,ct_4,ct_t];
Tip_x = h.*ones(1,3);

%% Geometric plot
plot(Root_x, Root_y,'b'); hold on;plot(Tip_x, Tip_y,'b');hold on;...
    axis equal; % cord plot
r_4  = [0 c_4; h ct_4]; % c/4 sweep plot
plot(r_4(:,1),r_4(:,2),'c'); hold on;

l_edge = [0 c; h ct_l]; % leading edge plot
plot(l_edge(:,1),l_edge(:,2),'k');hold on;

t_edge = [0 0; h ct_t]; % trailing edge plot
plot(t_edge(:,1), t_edge(:,2), 'k-'); hold on

%%   MAC geometry plot

I_rootl = c + t_c; I_roott = 0 - t_c; % plot imaginary tip cord (wingroot)
I_r1 = [0, I_rootl; 0, c]; plot(I_r1(:,1),I_r1(:,2),'g');hold on;
I_r2 = [0, 0; 0, I_roott]; plot(I_r2(:,1),I_r2(:,2),'g');hold on;

I_tipl = ct_l + c; I_tipt = ct_t - c; % plot imaginary root cord(wingtip)
I_t1 = [h, I_tipl; h, ct_l]; plot(I_t1(:,1),I_t1(:,2),'g'); hold on;
I_t2 = [h, ct_t; h, I_tipt]; plot(I_t2(:,1),I_t2(:,2),'g'); hold on;

% Intersecting lines
d_1 = [0, I_rootl; h, I_tipt]; plot(d_1(:,1),d_1(:,2),'g--'); hold on;
d_2 = [0, I_roott; h, I_tipl]; plot(d_2(:,1),d_2(:,2),'g--'); hold on;

%CG location
x = [0,h]; m = (-1.*(I_rootl-I_tipt))./h; y1 = (m.*x)+I_rootl;
m2 = (I_tipl-I_roott)./h; y2 = (m2.*x)+I_roott;...
    MAC_x = (I_roott-I_rootl)./(m-m2);...
    MAC_y = (m.*MAC_x)+I_rootl; plot(MAC_x,MAC_y,'ro');hold on;

% length MAC
MAC_l = c.*(2/3).*((1+t+(t.^2))./(1+t));
MAC_half = MAC_l/2;
MAC_lnew = [MAC_x, (MAC_y+MAC_half);MAC_x, (MAC_y-MAC_half)];
plot(MAC_lnew(:,1),MAC_lnew(:,2),'b'); hold on;

fprintf('Mean Aerodynamic Cord length (m): %i.\n',MAC_l)
fprintf('CG y-axis location: %i.\n',MAC_y)