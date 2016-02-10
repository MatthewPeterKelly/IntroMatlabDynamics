function MAIN_0_projectileMotion_MPK()

param.mass = 1.0;  % (kg)  mass of ball
param.gravity = 9.81;  % (m/s^2)  gravity accel.

% Initial state:
x0 = 0.0;  %initial position
y0 = 0.0;  %initial height
dx0 = 1.0; %initial x vel.
dy0 = 8.0;  %initial y vel.
z0 = [x0;y0;dx0;dy0];   % Initial state vector

% Set up for ODE45
time = linspace(0,2,100);
[time, state] = ode45(@ballDynamics,time,z0);

% Plot state vs time:
figure(1); clf;

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
figure(2); clf;
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
        
        % Equations of motion:  (compute accel.)
        ddx = 0;
        ddy = -g;    % -m*g = m*ddx
        
        % Pack up derivatives:
        dz = [ dx;dy;ddx;ddy];
        
    end




end

