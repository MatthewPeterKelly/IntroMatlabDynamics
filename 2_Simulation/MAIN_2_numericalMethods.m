function MAIN_2_numericalMethods()
%
% MAIN  --  basic numerical methods for integration
%
% Solve the simple harmonic oscillator equations using Euler's method, the
% mid-point method, Runge-Kutta, and ode45.
% 
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
%
% NOTES:
%
% I've used different time steps for all of the methods, with larger
% time-steps on the high-order methods. Play around with the time steps and
% see how the solutions behave.
%
% Euler's method cannot do a good job on this problem. It is said to be
% numerically unstable. Notice that the position and velocity estimates are
% getting farther and farther away from the true solution. A smaller
% time-step will help, but not fix the problem. This is because the
% integration is always adding a small amount of energy to the solution. 
%
% The mid-point method is only slightly harder to use than Euler's method,
% but it is a huge improvement! Notice that even with a much larger time
% step the solution does not diverge (Although if you make the time step
% really large then bad things will start to happen).
%
% "The" runge-kutta method is actually one of many Runge-Kutta methods. In
% fact, Euler's method is one of many first-order Runge-Kutta Methods. The
% Mid-point method is one of many second-order Runge-Kutta methods.
%
% Here we use the most commont 4th-order explicit runge-kutta method. It is
% a bit tricky to use and understand, but it is quite accurate, especially
% for systems with smooth dynamics. Notice that I've used a huge time step
% and it is still performing better than the mid-point method.
%
% Ode45 clearly does the best on the problem. Why? Well, let's start by
% looking at what ODE45 is doing. Notice that it has many time steps right
% at the beginning of the simulation. Any idea why? 
%
% It is because ode45 adaptively changes it's time step. Since it doesn't
% know what the scale of the problem is, it starts with a conservative
% estiamte (small step), and then computes an estimate of the error. If the
% estimate shows small error, then it will increase the time step, until it
% reaches the desired trade-off between cpu time and accuracy. This is what
% 'AbsTol' and 'RelTol' are controlling in the options struct.
%
% If ode45 is using variable grid points, then what happens when I pass in
% a uniform time vector? Answer: interpolating between it's own
% non-uniformly spaced grid points.
%
% Ode45 uses a high-order method, which is derived by cancelling error
% terms in a taylor series expansion. Look on Wikipedia if you're curious.
% It turns out that this method is the same as assuming that the simulation
% is a piece-wise high-order polynomial, and you can calculate this
% polynomial from the intermediate calculations that ode45 does. Try
% running ode45 with a single output (like done below). It will provide a
% solution struct, containing all of the data that describes this
% piece-wise polynomial solution, in sol.idata.f3d.
%
% When you pass in a uniform vector of times, it takes that piece-wise
% polynomial solution and evaluates it at the vector of times that you pass
% in.
%
%



% Initial state (position and velocity)
x0 = 1.0;   % Initial position
v0 = 0.0;   % Initial velocity
z0 = [x0;v0];   % Initial state

% Time span for the simulation:
tSpan = [0,15];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         Analytic Solution                               %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Define a function handle for the position solution
tSoln = linspace(tSpan(1), tSpan(2), 100);
xSoln = v0*sin(tSoln) + x0*cos(tSoln);
vSoln = v0*cos(tSoln) - x0*sin(tSoln);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                          Euler's Method                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Define the solution grid for Euler's method solution
dtEuler = 0.06;
tEuler = tSpan(1) : dtEuler : tSpan(2);

% Call numerical integration with Euler's method:
zEuler = eulerMethod(@simpleHarmonicOscillator,tEuler,z0);
xEuler = zEuler(1,:);
vEuler = zEuler(2,:);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Midpoint Method                                  %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Define the solution grid for Midpoint method solution
dtMidPt = 0.3;
tMidPt = tSpan(1) : dtMidPt : tSpan(2);

% Call numerical integration with MidPoint method:
zMidPt = midPtMethod(@simpleHarmonicOscillator,tMidPt,z0);
xMidPt = zMidPt(1,:);
vMidPt = zMidPt(2,:);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                    "The" Runge-Kutta Method                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Define the solution grid for Runge-Kutta method solution
dtRK4 = 1;
tRK4 = tSpan(1) : dtRK4 : tSpan(2);

% Call numerical integration with MidPoint method:
zRK4 = rungeKutta(@simpleHarmonicOscillator,tRK4,z0);
xRK4 = zRK4(1,:);
vRK4 = zRK4(2,:);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                          Ode45 Method                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Call numerical integration with ode45 method. Return the full solution
% struct, rather than the interpolant of the solution
soln = ode45(@simpleHarmonicOscillator,tSpan,z0);

% Solution on the time-stepping grid that ode45 actually used
tOde45Grid = soln.x;   %Time (which ode45 confusingly calls 'x')
zOde45Grid = soln.y;   %State (which ode45 confusingly calls 'y')
xOde45Grid = zOde45Grid(1,:);  % Extract position
vOde45Grid = zOde45Grid(2,:);  % Extract velocity

% Interpolate ode45 solution
tOde45 = linspace(tSpan(1), tSpan(2), 100);
zOde45 = deval(soln,tOde45);    %Fancy method-consistent interpolation
xOde45 = zOde45(1,:);
vOde45 = zOde45(2,:);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Compare Solutions                                %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

figColor = 0.7*[1,1,1];
ode45Color = [0.2,0.8,0.2];
eulerColor = [0.8,0.2,0.2];
midPtColor = [0.2,0.2,0.8];
rk4Color = [0.6,0.2,0.7];

figure(10); clf; set(gcf,'color',figColor)

subplot(2,1,1); hold on; set(gca,'color',figColor)
plot(tSoln,xSoln,'k-','LineWidth',3)
plot(tEuler,xEuler,'o',...
    'LineWidth',1,'MarkerSize',4,'Color',eulerColor);
plot(tMidPt,xMidPt,'o',...
    'LineWidth',1,'MarkerSize',6,'Color',midPtColor);
plot(tRK4,xRK4,'o',...
    'LineWidth',2,'MarkerSize',7,'Color',rk4Color);
plot(tOde45,xOde45,'-','LineWidth',1,'Color',ode45Color)
plot(tOde45Grid,xOde45Grid,'o',...
    'LineWidth',3,'MarkerSize',10,'Color',ode45Color);

xlabel('time (s)')
ylabel('position (m)')
title('position vs time')
legend('analytic','euler','midpoint','rk4','ode45','ode45 grid',...
    'Location','NorthEastOutside');

subplot(2,1,2); hold on; set(gca,'color',figColor)
plot(tSoln,vSoln,'k-','LineWidth',3)
plot(tEuler,vEuler,'o',...
    'LineWidth',1,'MarkerSize',4,'Color',eulerColor);
plot(tMidPt,vMidPt,'o',...
    'LineWidth',1,'MarkerSize',6,'Color',midPtColor);
plot(tRK4,vRK4,'o',...
    'LineWidth',2,'MarkerSize',7,'Color',rk4Color);
plot(tOde45,vOde45,'-','LineWidth',1,'Color',ode45Color)
plot(tOde45Grid,vOde45Grid,'o',...
    'LineWidth',3,'MarkerSize',10,'Color',ode45Color);

xlabel('time (s)')
ylabel('velocity (m)')
title('velocity vs time')
legend('analytic','euler','midpoint','rk4','ode45','ode45 grid',...
    'Location','NorthEastOutside');
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


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


function z = eulerMethod(dynFun,t,z0)
% z = eulerMethod(dynFun,t,z0)
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



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


function z = midPtMethod(dynFun,t,z0)
% z = midPtMethod(dynFun,t,z0)
% 
% Simulates a system using midPtMethod's method
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
    
    % Compute the state and dynamics at lower edge of the interval
    tLow = t(i-1);   %Current time
    zLow = z(:,i-1);  %Current state
    dzLow = dynFun(tLow,zLow);  %Derivative at this point:
   
    % Estimate the state and dynamics at the mid-point:
    tMid = tLow + 0.5*dt;
    zMid = zLow + 0.5*dt*dzLow;
    dzMid = dynFun(tMid,zMid);
    
    % Use the mid-point dynamics estimate to compute simulation step:
    zUpp = zLow + dt*dzMid;
    z(:,i) = zUpp; 
end

end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


function z = rungeKutta(dynFun,t,z0)
% z = rungeKutta(dynFun,t,z0)
% 
% Simulates a system using "The Runge-Kutta" method (4th-order)
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
% NOTES:
%
%   This is a bit like the mid-point method, but with a few more
%   intermediate steps. This is a 4th-order method, and thus is much more
%   accurate than a 2nd-order method like the Mid-point method, or a
%   1st-order method like Euler's method.
%
%   ODE45 uses a very similar method. As the name implies, ode45 uses a
%   combined 4th and 5th order method. Since it is running two methods at
%   once, it can compute an estimate of the error at the end of the step by
%   comparing the 5th-order method with the 4th-order method. It uses this
%   error estimate to adjust the step size on the fly. Notice that all of
%   the methods that I've written out here are "fixed-step" methods: the
%   just march along a uniform grid. ode45 does not. It uses a variable
%   grid, that it computes on the fly. When you pass a uniform time grid to
%   ode45, it is actually using high-order interpolation between the
%   non-uniform grid points of the solution. This interpolation scheme is
%   derived along with the integration scheme, and is thus also 4th-order
%   accurate.
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
    
    % Compute the dynamics estimate at the lower edge of the interval    
    tLow = t(i-1);   %Current time
    zLow = z(:,i-1);  %Current state
    k1 = dynFun(tLow,zLow);  
    
    % Use the previous dynamics estimate to march forward to the mid-point
    tMid = t(i-1) + 0.5*dt;
    zMidA = zLow + 0.5*dt*k1;
    k2 = dynFun(tMid, zMidA);
    
    % Now compute a better estimate of the mid-point state and dynamics
    zMidB = zLow + 0.5*dt*k2;
    k3 = dynFun(tMid,zMidB);
    
    % Finally, estimate the upper edge of the segment:
    tUpp = tLow + dt;
    zUpp = zLow + dt*k3;
    k4 = dynFun(tUpp,zUpp);
    
    % This part is a bit tricky. Estimate the integral using a weighted sum
    % of the dynamics at each of the estimated points. This equation is not
    % arbitrary, and it is calculated to cancel low-order error terms.
    zUpp = zLow + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    z(:,i) = zUpp; 
end

end

