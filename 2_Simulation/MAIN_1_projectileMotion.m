function MAIN_1_projectileMotion()
%
% MAIN  --  Projectile motion using ode45
%
% 2D point mass in constant gravity field:
% Basic simulation in Matlab using ode45 with 2 D.O.F
%

% initialize values
x0 = 0; y0 = 10; % initial x & y position of particle [m]
vx0 = 5; vy0 = 10; % initial x & y velocities of particle [m/s]
z0= [x0,y0,vx0,vy0]'; % Pack ICs into z0 and make it a column vector
t = linspace(0,2,80); % times to look at solution [sec] (0 is initial
% time, 2 is final time, with 100 evenly spaced
% points)

% Parameters are stored in 'p' structure for easy bookkeeping
p.m = 1; % set point mass to be 1 [kg]
p.g = 9.81; % set acceleration due to gravity [m/s^2]

% ODESET sets preferences for ODE solvers like ODE45
% here we tell ode45 how much error we tolerate in the solution
options = odeset('abstol', 1e-6,'reltol',1e-6);

% The next line numerically solves the ODEs that are in massInGravityField
[~,z] = ode45(@massInGravityField,t, z0,options,p);

% Pull out x and y comps of position at different times
x = z(:,1);
y = z(:,2);

% Plot the trajectory (x vs y)
figure(1); clf;
hold on; % Retain current plot when performing following commands
axis equal; % Use the same length for the data units along each axis
plot(x,y,'k--','linewidth',2)
xlabel('x-position (m)','fontsize',14)
ylabel('y-position (m)','fontsize',14)
title('Trajectory of Point Mass','fontsize',16)

% Animation
for i=1: length(t);
    if i==1
        % Create a marker for the projectile and plot it on the figure
        circleHandle = plot(x(), y(i),...   
            'ro','markersize',14,'linewidth',6);
    else
        % Now that marker is created, just move it on each new frame 
        set(circleHandle,'xData',x(i),'yData',y(i));
    end
    pause(.05);    %Pause between frames to make animation visible
    drawnow;   % Force matlab to draw the figure
end

% Plot the solution against time:
figure(2); clf;

subplot(2,1,1);
plot(t,x)
xlabel('time (s)')
ylabel('x-position (m)')
title('distance vs time')

subplot(2,1,2);
plot(t,y)
xlabel('time (s)')
ylabel('y-position (m)')
title('height vs time')


end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function dz = massInGravityField(t,z,p)
%
% Computes the dynamics a particle in a constant gravitational field
%
% INPUTS:
%   t = current time
%   z = [x;y;dx;dy] = state vector
%   p = structure containing parameters
%
% OUTPUTS:
%   dz = [dx;dy;ddx;ddy] = time-derivative of the state vector
%

% unpack parameters
m = p.m; % particle mass [kg]
g = p.g; % acceleration due to gravity [m/s^2]

x=z(1); y= z(2); dx = z(3);  dy = z(4);  % unpack z into readable names

% x and y components of force on particle
Fx = 0; % no horizontal force component
Fy = -m*g; % vertical force component from gravity

% compute x and y components of acceleration using Newton's 2nd (a=F/m)
ddx = Fx/m; % x acceleration component
ddy = Fy/m; % y acceleration component

% The rates of change of all variables in z
dz = [dx,dy,ddx,ddy]';
end