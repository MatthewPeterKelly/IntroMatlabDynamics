function MAIN_0_Euler()
%
% MAIN  --  Euler
%
% Basic simulation in Matlab, using vectors, functions, and Euler's method
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

% Time step for the numerical integration method. Notice that the solution
% by Euler's method is not very good when the time step is large!
h = 0.05;   %Time step 

% Set up timing stuff
t0 = 0;
tF = 10;
t = t0 : h : tF;  %Time vector



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         Analytic Solution                               %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


% Initial state (position and velocity)
x0 = 1.0;
v0 = 0.0;

% Define a function handle for the solution:
xSoln = @(t)( v0*sin(t) + x0*cos(t) );
vSoln = @(t)( v0*cos(t) - x0*sin(t) );

% Make a function for computing the state vector:
zSoln = @(t)( [xSoln(t); vSoln(t)] );

% Plot the solution:
figure(1); clf;

subplot(2,1,1);
plot(t,xSoln(t),'k-');
xlabel('time (s)')
ylabel('position (m)')
title('Position vs Time')

subplot(2,1,2);
plot(t,vSoln(t),'k-');
xlabel('time (s)')
ylabel('velocity (m/s)')
title('Velocity vs Time')


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%             Numerical Solution      (Euler's Method)                    %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Initial State vector
z0 = [x0;v0];

% Call numerical integration with Euler's method:
dynFun = @(t,z)( simpleHarmonicOscillator(z) );
z = eulerMethod(dynFun,t,z0);
x = z(1,:);
v = z(2,:);

% Plot the numerical solution over top of the analytic solution:
subplot(2,1,1); hold on
plot(t,x,'ro');

subplot(2,1,2); hold on;
plot(t,v,'ro');


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                Numerical Solution    (ode45)                            %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Initial State vector
z0 = [x0;v0];

% Call numerical integration with Euler's method:
dynFun = @(t,z)( simpleHarmonicOscillator(z) );
[~, z] = ode45(dynFun,t,z0); z = z';
x = z(1,:);
v = z(2,:);

% Plot the numerical solution over top of the analytic solution:
subplot(2,1,1); hold on
plot(t,x,'bs');
legend('analytic','euler','ode5')

subplot(2,1,2); hold on;
plot(t,v,'bs');
legend('analytic','euler','ode5')


end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function dz = simpleHarmonicOscillator(z)
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



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function z = eulerMethod(dynFun,t,z0)
%
% Simulates a system using Euler's method
%
% INPUTS:
%   dynFun = function handle
%       dz = dynFun(t,z)
%           t = [1,1] = time
%           z = [nz,1] = state vector
%           dz = [nz,1] = dz/dt = derivative of state vector
%   t = [1,nt] = time vector
%   z0 = [nz,1] = initial state
% 
% OUTPUTS:
%   z = [nz, nt] = state at each grid point
%

% Figure out problem size
nt = length(t);
nz = size(z0,1);

% Allocate memory for the output:
z = zeros(nz,nt);

% Store the initial state:
z(:,1) = z0;

% March forward in time:
for i=2:nt
    dt = t(i)-t(i-1);  %Current time step
    tNow = t(i-1);   %Current time
    zNow = z(:,i-1);  %Current state
   
    dzNow = dynFun(tNow,zNow);  %Derivative at this point:
   
    z(:,i) = zNow + dt*dzNow; 
end

end









