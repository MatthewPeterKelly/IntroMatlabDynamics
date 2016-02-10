function MAIN_1_twoSimplePendulums()
%
% MAIN  --  An ode45 simulation of 2 separate simple pendulums of different
% lengths.
%
% Equations Of Motion:
%  thetadotdot = -(g/L)sin(theta)
%
% theta1 = angle of the first pendulum
% omega1 = theta1dot = angular velocity of first pendulum
% omega1dot = angular acceleration of first pendulum
% theta2 = angle of the second pendulum
% omega2 = theta2dot = angular velocity of second pendulum
% omega2dot = angular acceleration of second pendulum
%

% initialize values
theta10 = pi-0.1; theta20 = pi-0.1; % initial angles of pendulums [rad]
omega10 = 0; omega20 = 0; % initial angular velocities of pendulums [rad/s]

% Pack ICs into initial state vector z0 and make it a column vector
z0= [theta10,theta20,omega10,omega20]'; 

time = linspace(0,10,200); % times to look at solution [sec] (0 is initial
% time, 10 is final time, with 200 evenly spaced points)

% Parameters are stored in 'p' structure for easy bookkeeping
p.L1 = 1; p.L2 = 1.5; % set lengths of pendulum arms [m]
p.g = 9.81; % set acceleration due to gravity [m/s^2]

% Create a nested function for the dynamics. This allows us to pass the
% parameter struct (p) as a special type of global variable, which avoids
% the confusing (and undocumented) syntax for ode45 where you tack on
% parameters to the end of the argument list in ode45.
    function zdot = twoSimplePendulums(t,z)
        %
        % Computes the dynamics a particle in a constant gravitational field
        %
        % INPUTS:
        %   t = current time
        %   z = [theta1;theta2;omega1;omega2] = state vector
        %
        % OUTPUTS:
        %   zdot = [theta1dot;theta2dot;omega1dot;omega2dot] = time-derivative of the state vector
        %
        
        % unpack parameters
        L1 = p.L1; L2 = p.L2; % particle mass [kg]
        g = p.g; % acceleration due to gravity [m/s^2]
        
        % unpack z into readable names
        theta1=z(1); theta2= z(2);
        omega1 = z(3);  omega2 = z(4);  
        
        % Compute time-derivative of the state vector using equations of motion:
        theta1dot = omega1;
        theta2dot = omega2;
        omega1dot = -(g/L1)*sin(theta1);
        omega2dot = -(g/L2)*sin(theta2);
        
        % The rates of change of all variables in z
        zdot = [theta1dot;theta2dot;omega1dot;omega2dot];
    end

% ODESET sets preferences for ODE solvers like ODE45
% here we tell ode45 how much error we tolerate in the solution
options = odeset('abstol', 1e-6,'reltol',1e-6);

% The next line numerically solves the ODEs that are in twoSimplePendulums
[~,zarray] = ode45(@twoSimplePendulums,time,z0,options);

% Pull out theta1 and theta2 comps of position at different times
theta1soln = zarray(:,1);
theta2soln = zarray(:,2);

% Plot the solution against time:
figure(1); clf;

subplot(2,1,1);
plot(time,theta1soln)
xlabel('time (s)')
ylabel('theta1 (rad)')
title('angle vs time')

subplot(2,1,2);
plot(time,theta2soln)
xlabel('time (s)')
ylabel('theta2 (rad)')
title('angle vs time')

% Animation of the first pendulum
x1=L1*sin(theta1soln); %convert polar coordinates to cartesian
y1=-L1*cos(theta1soln);
%create a figure for the animation
figure(2); clf;
hold on % plot everything on the same plot
axis equal % units same length along each axis
% set axis limits to contain entire animation
axis([min(x1),max(x1),min(y1),max(y1)]);
for i=1: length(time);
    if i==1
        % Create a marker for the point mass and plot it on the figure
        circleHandle = plot(x1(i),y1(i),...
            'ko','markersize',14,'linewidth',6);
        % Create line representing the massless string
        lineHandle = plot([0,x1(i)],[0,y1(i)],...
            'k','linewidth',3);
    else
        % Now that marker and line are created, just move it on each new frame
        set(circleHandle,'xData',x1(i),'yData',y1(i));
        set(lineHandle,'xData',[0,x1(i)],'yData',[0,y1(i)]);
    end
    pause(.05);    %Pause between frames to make animation visible
    drawnow;   % Force matlab to draw the figure
end

end
