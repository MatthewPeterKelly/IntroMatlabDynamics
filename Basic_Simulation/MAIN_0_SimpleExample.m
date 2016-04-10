function MAIN_0_SimpleExample()
%
% Simple example for ode45
%
% System:  ddx = x + dx
%
%
% There are a few basic steps to simulation a system in Matlab.
%
% 1) Write down the system dynamics as a function that can then be passed
% into ode45 using the @ symbol to create a function handle
%
% 2) Set up the simulation. You need to tell ode45 what the initial state
% of your system is, when it should start and stop the simulation, and what
% times you would like to evaluate the solution at. You can also pass in an
% options struct for more details (>> help ode45)
%
% 3) Once you have the solution from ode45, you need to unpack it to get the
% original coordinates out of your state vector. In this case, we want to
% get the position and velocity.
%
% 4) The last step is to do something with the solution, usually this
% involves either plotting or animating the system.
%

% Start by constructing a nested function that describes your dynamics:
    function dz = dynamicsFunction(t,z)
        % Computes first-order form of the system dynamics:
        % ddx = x + dx
        %
        % INPUTS:
        %   t = scalar time = [not used]
        %   z = [2,1] = [x;v] = state (column) vector
        %   
        % OUTPUS
        %   dz = dz/dt = time derivative of state
        %
        
        % Unpack the state
        x = z(1);
        v = z(2);
        
        % Compute (first-order form) dynamics:
        dx = v;
        dv = x + v;
        
        % Pack back up into derivative vector:
        dz = [dx;dv];
        
    end

% next, set up the simulation time span and time vector
tSpan = [0,2];  
tSoln = linspace(tSpan(1), tSpan(2), 100);

% Initial state of the system:
x0 = 0.1;
v0 = -0.2;
z0 = [x0;v0];

%%%% now call ode45:
% 1) The first argument is a function handle, which is created using @
% 2) The second argument is a time span, or vector of times at which you
% would like the state of simulation to be evaluated
% 3) The third argument is the initial state of the system
[~, z] = ode45(@dynamicsFunction,tSoln,z0);  z = z';

% Unpack the solution
xSoln = z(1,:);
vSoln = z(2,:);

% Plots:
figure(1); clf;

subplot(2,1,1);
plot(tSoln, xSoln)
title('position vs time')
xlabel('time (s)')
ylabel('position (m)')

subplot(2,1,2);
plot(tSoln, vSoln)
title('velocity vs time')
xlabel('time (s)')
ylabel('velocity (m/s)')

end