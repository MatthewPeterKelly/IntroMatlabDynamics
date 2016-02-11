function MAIN_2a_projectileMotionWithDrag()
%
% This function uses ode45 to simulate a ball flying through the air,
% subject to linear and quadratic drag forces.
%
% ![Advanced Version]!
%
% Shows how to use a function handle to pass function parameters to a
% sub-function, instead of using a nested function. For details, try:
%   >> help function_handle
%
% This also shows how to check which points ode45 actually uses, and then
% how to interpolate the solution between them.
%

% Struct with problem parameters:
param.mass = 1.0;  % (kg)  mass of ball
param.gravity = 9.81;  % (m/s^2)  gravity accel.
param.linearDragCoefficient = 0.8;  % (Ns/m)
param.quadraticDrawCoefficient = 0.7; % (Ns^2/m^2)

% Initial state:
x0 = 0.0;  %initial position
y0 = 0.0;  %initial height
dx0 = 4.0; %initial x vel.
dy0 = 6.0;  %initial y vel.
z0 = [x0;y0;dx0;dy0];   % Initial state vector

% Set up for ODE45
tSpan = [0, 0.8];
dynFun = @(t,z)( ballDynamics(z,param) );  %Anonymous function (lambda function)
soln = ode45(dynFun,tSpan,z0);  %Return the solution struct

% unpack the solution on the actual grid-points used by ode45:
tGrid = soln.x;
zGrid = soln.y;
xGrid = zGrid(1,:);
yGrid = zGrid(2,:);

% Interpolate the solution on a uniform grid, for plotting. This is what
% ode45 normally does inside, when you pass in a time vector.
t = linspace(tSpan(1), tSpan(2), 100);
z = deval(soln, t);  %high-order accurate interpolation for ode solvers
x = z(1,:);
y = z(2,:);

% Plot state vs time:
figure(4); clf;

subplot(2,1,1); hold on;
plot(tGrid,xGrid,'ko','MarkerSize',8,'LineWidth',1);
plot(t,x);
xlabel('time (s)')
ylabel('x position (m)');
title('Projectile Motion!')
legend('ode45 grid', 'interpolant');

subplot(2,1,2); hold on;
plot(tGrid,yGrid,'ko','MarkerSize',8,'LineWidth',1);
plot(t,y);
xlabel('time (s)')
ylabel('y position (m)');
legend('ode45 grid', 'interpolant');

% Plot trajectory shape
figure(5); clf; 
hold on;
plot(xGrid,yGrid,'ko','MarkerSize',8,'LineWidth',1);
plot(x,y);
xlabel('x position (m)');
ylabel('y position (m)');
title('Trajectory Shape');
legend('ode45 grid', 'interpolant');
axis equal;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dz = ballDynamics(z, param)
%
% Computes the dynamics for a ball flying through the air with linear and
% quadraticc drag.
%
% INPUTS:
%   z = [x;y;dx;dy] = state vector
%
% OUTPUTS: 
%   dz = dz/dt = time derivative of the state
%

% Unpack state:
x = z(1); % Horizontal position of ball
y = z(2); % Vertical position of ball
dx = z(3); % horiz. velocity
dy = z(4); % vert. velocity

% Unpack parameters:
m = param.mass;
g = param.gravity;
c1 = param.linearDragCoefficient;
c2 = param.quadraticDrawCoefficient;

% Drag force:  (Linear + Quadratic)
v = sqrt(x^2 + y^2);   % speed of the ball
Fdx = (-c1*dx) + (-c2*dx*v);   % horizontal terms
Fdy = (-c1*dx) + (-c2*dx*v);   % vertical terms

% Equations of motion:
ddx = Fdx/m;            %vertical acceleration
ddy = (Fdy - m*g)/m;    %horizontal acceleration

% Pack up derivatives:
dz = [dx; dy; ddx; ddy];

end
