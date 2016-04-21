%%%% Simulation of the disk rolling on the see-saw
%
%

% Initial State
z0 = [...
    0.1;  % absolute angle of the see-saw
    0.6;  %absolute angle of the disk
    0.0;  % absolute angle rate of the see-saw
    -2.6]; % absolute angle rate of the disk

% Physical Parameters
p.m = 1.0;  %(kg)mass of the disk
p.R = 0.2;  %(m) radius of the disk
p.I2 = 0.5*p.m*p.R^2;  %moment of inertia of the disk about its CoM
p.I1 =  10*p.I2;  %inertia of the see-saw about CoM == pivot
p.g = 9.81; % (m/s^2);   %gravity acceleration

% Set up and solve simulation:
tSpan = [0,2];  
dynFun = @(t,z)( dynamics(z,p) );
options = odeset(...
    'AbsTol',1e-12,...
    'RelTol',1e-12);
sol = ode45(dynFun,tSpan,z0,options);

% Unpack the solution:
t = linspace(tSpan(1), tSpan(2), 100);
z = deval(sol,t);
q1 = z(1,:);
q2 = z(2,:);
dq1 = z(3,:);
dq2 = z(4,:);
[E,T,U] = energy(z,p);

% Plots:
figure(1); clf;

subplot(2,3,1);
plot(t,q1);
xlabel('time');
ylabel('angle');
title('See-Saw Angle')

subplot(2,3,2);
plot(t,q2);
xlabel('time');
ylabel('angle');
title('Disk Angle')

subplot(2,3,3); hold on;
plot(t,E)
plot(t,T)
plot(t,U)
xlabel('time')
ylabel('energy')
title('Energy');
legend('total','kinetic','potential');

subplot(2,3,4);
plot(t,dq1);
xlabel('time');
ylabel('rate');

subplot(2,3,5);
plot(t,dq2);
xlabel('time');
ylabel('rate');
