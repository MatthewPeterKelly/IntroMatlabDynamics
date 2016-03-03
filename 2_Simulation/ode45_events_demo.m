function [t,y,dy] = ode45_events_demo()
%
% Show how to use ode45 to compute a simulation with events, in this case:
%
% A ball is thrown straight up. Compute a simulation of it's trajectory, in
% one dimension, given that it is influenced by gravity and quadratic air
% drag. Stop the simulation exactly when the ball hits the ground.
%

p.g = 1;    %Acceleration due to gravity
p.c = 0.1;   % quadratic drag coefficient

y0 = 5;  % Initial ball height
dy0 = 2;   % initial ball velocity

z0 = [y0; dy0];  % Initial state
tSpan = [0,10];  % [ initial time,  maximum time ] 

options = odeset(...   % dots --> tells matlab to read onto the next line
    'Events',@(t,z)( groundEvent(z) )...
    );  %Function handle for events

dynFun = @(t,z)( dynamics(z,p) );   % Pull struct p from current workspace

soln = ode45(dynFun,tSpan,z0,options);  % soln = solution struct

collisionTime = soln.xe;  % Time at which event occured
t = linspace(0, collisionTime, 100);  % time vector for plotting
z = deval(soln,t);  % Interpolate solution at desired times
y = z(1,:);  %height
dy = z(2,:); %velocity

%Plotting:
figure(10); clf;
subplot(2,1,1);
plot(t,y)
xlabel('time')
ylabel('height')
title('Ball falling to ground')
subplot(2,1,2);
plot(t,dy)
xlabel('time')
ylabel('speed')

% Print out:
fprintf('Collision occured at time = %4.4f \n',collisionTime);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dz = dynamics(z,p)
%
% Computes the dynamics of a point mass falling under influence of gravity
% and quadratic drag.
%
% INPUTS:
%   z = [y;dz] = [height; speed];
%   p = struct of parameters
%       .g = gravity acceleration
%       .c = drag coeff
%
% OUTPUTS:
%   dz = dz/dt = [speed; acceleration];
%

% y = z(1,:);   %Not needed
dy = z(2,:);

% acceleration = gravity + drag
ddy = (-p.g) + (-p.c*dy);

% Collect terms:
dz = [dy;ddy];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value, isTerminal, direction] = groundEvent(z)

y = z(1,:);   
% dy = z(2,:);  %Not needed

value = y;   %Event occurs when value == 0  (ball hits ground)
isTerminal = true;  %ode45 should stop when event occurs
direction = -1; %Sign of derivative of value at event  (0 would also work)

%%%% A note about direction:
% direction = 1  --> Only trigger event if derivative of value is positive
% direction = -1 --> Only trigger event if derivative of value is negative
% direction = 0  --> Always trigger event
%
% why?
% --> Event detection is computationally expensive, and there are some
% applications where you only care about events with a specific derivative.
%

end


