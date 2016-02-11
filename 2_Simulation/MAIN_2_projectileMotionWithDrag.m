function MAIN_2_projectileMotionWithDrag()
%
% This function uses ode45 to simulate a ball flying through the air,
% subject to linear and quadratic drag forces. 
%

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
time = linspace(0, 0.8, 100);
[time, state] = ode45(@ballDynamics,time,z0);

% Plot state vs time:
figure(4); clf;

subplot(2,1,1);
plot(time,state(:,1))
xlabel('time (s)')
ylabel('x position (m)');
title('Projectile Motion!')

subplot(2,1,2); 
plot(time,state(:,2))
xlabel('time (s)')
ylabel('y position (m)');


% Plot trajectory:
figure(5); clf;
xlabel('x position (m)');
ylabel('y position (m)');
title('Trajectory Shape');

hold on;
plot(state(:,1),state(:,2));
plot(state(1,1),state(1,2),'rx',...
    'LineWidth',4,'MarkerSize',15);
axis equal;


    function dz = ballDynamics(t,z)
                
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
        dz = [ dx;dy;ddx;ddy];
        
    end


end

