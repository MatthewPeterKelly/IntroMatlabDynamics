function MAIN_0_simpleHarmonicOscillator()
%
% MAIN  --  Use ODE45 to simulate a simple harmonic oscillator
%
% Basic simulation in Matlab using ode45
%
% Simple Harmonic Oscillator:
% 
%   t = time (s)
%   x = position (m)
%   dx = velocity (m/s)
%   ddx = acceleration (m/s^2)
%
%   Equations of motion:
%   ddx + x = 0
%
%   First order form:   (set v = dx);
%   dx = v;
%   dv = -x;
%
%   Solution:
%   x(0) = x0;
%   v(0) = v0;
%   x(t) = v0*sin(t) + x0*cos(t)
%   
%   State vector:
%   z = [x;v]
%

% Time span for the simulation
tSpan = [0,1];

% Initial state (position and velocity)
x0 = 1.0;
v0 = 0.0;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         Analytic Solution                               %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


% Analytic Solution:
tSoln = linspace(tSpan(1), tSpan(2), 100);
xSoln = v0*sin(tSoln) + x0*cos(tSoln);
vSoln = v0*cos(tSoln) - x0*sin(tSoln);

% Plot the solution:
figure(1); clf;

subplot(2,1,1);
plot(tSoln,xSoln,'k-');
xlabel('time (s)')
ylabel('position (m)')
title('Position vs Time')

subplot(2,1,2);
plot(tSoln,vSoln,'k-');
xlabel('time (s)')
ylabel('velocity (m/s)')
title('Velocity vs Time')


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                Numerical Solution    (ode45)                            %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Vector of times that we would like ode45 to return the state at:
t = linspace(tSpan(1), tSpan(2), 100);

% Initial State vector
z0 = [x0;v0];

% Call numerical integration with Euler's method:
[~, z] = ode45(@simpleHarmonicOscillator,t,z0); z = z';
x = z(1,:);
v = z(2,:);

% Plot the numerical solution over top of the analytic solution:
subplot(2,1,1); hold on
plot(t,x,'bs');
legend('analytic','ode45')

subplot(2,1,2); hold on;
plot(t,v,'bs');
legend('analytic','ode5')


end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function dz = simpleHarmonicOscillator(t,z)
%
% Computes the dynamics of the simple harmonic oscillator
%
% INPUTS:
%   z = [2, n] = [x;v] = state vector
%       z(1,:) = x = position
%       z(2,:) = v = velocity
%
% OUTPUTS:
%   dz = dz/dt = time-derivative of the state vector
%

% Unpack State
x = z(1,:);
v = z(2,:);

% Compute dynamics
dx = v;
dv = -x;

% Pack up dynamics
dz = [dx;dv];

end







