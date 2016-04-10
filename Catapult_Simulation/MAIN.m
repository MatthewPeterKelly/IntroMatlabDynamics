% MAIN  --  Catapult Simulation
%
% Simulated a simple torsion catapult launching a projectile.
%
% Uses ode45 for simulation and event-detection.
%
% Phase One:  simulation of the catapult launching the projectil
% Phase Two:  simulation of the projectile flying through the air
%
clc; clear;
clear drawCatapultSimulation.m  % clears persistent variables

%%%% Problem parameters
param.armMass = 10;   %(kg)  arm is a slendar rod
param.projectileMass = 2;  %(kg)
param.armLength = 3;  %(m)
param.gravity = 9.81;  %(m/s^2)
param.springConstant = 1500;  %(N/rad)
param.springRestAngle = 0*(pi/180);  % (rad)  measured from pos. vert. axis.
param.initialAngle = (90+30)*(pi/180);  % (rad)  measured from pos. vert. axis.
param.quadraticAirDrag = 0.1;  %(N-s^2/m^2)
param.launchAngle = 45*(pi/180);   %(rad) measured from pos. vert. axis.
param.xCatapult = 0;  %(m)  horizontal position of catapult axle
param.yCatapult = 2;  %(m)  height of the catapult axle above ground



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                  ode45 simulation of launch                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
z0Launch = [param.initialAngle; 0];
tSpanLaunch = [0,4];
dynLaunch = @(t,z)( catapultDynamics(z,param) );
odeOptLaunch = odeset(...
    'Event',@(t,z)( launchEvent(z,param) ),...
    'AbsTol',1e-8,...
    'RelTol',1e-8);
solLaunch = ode45(dynLaunch,tSpanLaunch,z0Launch,odeOptLaunch);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                analysis of the launch simulation                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%% Check to make sure that the catapult launched
zLaunchEvent = solLaunch.ye;   % state of the system at the instant of launch
if isempty(zLaunchEvent)
    error('No launch event detected. Invalid parameter set.');
end

%%%% Extract the grid points that ode45 actually used:
tLaunchGrid = solLaunch.x;
qLaunchGrid = solLaunch.y(1,:);
dqLaunchGrid = solLaunch.y(2,:);

%%%% Interpolat solution between grid-points
tLaunch = linspace(tLaunchGrid(1), tLaunchGrid(end), 50);  %time grid for plotting results
zLaunch = deval(solLaunch,tLaunch);  % state at each grid point in time
qLaunch = zLaunch(1,:);  % angle of the catapult arm
dqLaunch = zLaunch(2,:); % angular rate of the catapult arm


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                  ode45 simulation of flight                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
z0Flight = getProjectileState(zLaunch(:,end), param);
collisionDuration = 1e-6;
tSpanFlight = tLaunch(end) + collisionDuration + [0,60];
dynFlight = @(t,z)( projectileDynamics(z,param) );
odeOptFlight = odeset(...
    'Event',@(t,z)( groundEvent(z) ),...
    'AbsTol',1e-8,...
    'RelTol',1e-8);
solFlight = ode45(dynFlight,tSpanFlight,z0Flight,odeOptFlight);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                analysis of the flight simulation                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%% Check to make sure that the catapult launched
zGroundEvent = solFlight.ye;   % state of the system at the instant of launch
if isempty(zGroundEvent)
    error('No ground event detected. Invalid parameter set.');
end

%%%% Extract the grid points that ode45 actually used:
tFlightGrid = solFlight.x;
xFlightGrid = solFlight.y(1,:);
yFlightGrid = solFlight.y(2,:);
dxFlightGrid = solFlight.y(3,:);
dyFlightGrid = solFlight.y(4,:);

%%%% Interpolat solution between grid-points
tFlight = linspace(tFlightGrid(1), tFlightGrid(end), 150);  %time grid for plotting results
zFlight = deval(solFlight,tFlight);  % state at each grid point in time
xFlight = zFlight(1,:);
yFlight = zFlight(2,:);
dxFlight = zFlight(3,:);
dyFlight = zFlight(4,:);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%               Compute full-state across simulation                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% zFull = [x;y;q;dx;dy;dq]

%%%% Compute projectile state during launch:
zTmp = getProjectileState(zLaunch, param);
xLaunch = zTmp(1,:);
yLaunch = zTmp(2,:);
dxLaunch = zTmp(3,:);
dyLaunch = zTmp(4,:);

%%%% Catapult state during flight  (constant, at hard stop)
qFlight = qLaunch(end)*ones(size(tFlight));
dqFlight = zeros(size(tFlight));

%%%% Pack everything up:
tFull = [tLaunch, tFlight];
qFull = [qLaunch, qFlight];
xFull = [xLaunch, xFlight];
yFull = [yLaunch, yFlight];
dxFull = [dqLaunch, dqFlight];
dyFull = [dxLaunch, dxFlight];
dqFull = [dyLaunch, dyFlight];
zFull = [xFull;yFull;qFull;dxFull;dyFull;dqFull];

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%               Animation of the catapult simulation                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

figure(4); clf; hold on;

% Colors and markers
lineWidth = 2;
lineColor = [0.2, 0.2, 0.8];  % Blue
markerColor = [0.8,0.2,0.2];  % Red
groundWidth = 5;
groundColor = [77,38,0]/255;   % Brown

% Get the scaling for the graphics
xBnd = [min(xFull), max(xFull)] + 0.2*[-1,1]*range(xFull);
xGround = linspace(xBnd(1), xBnd(2), 150);
yGround = groundModel(xGround);
yAll = [yFull,yGround];  %Collect all y points, used for scaling only
yBnd =  [min(yAll), max(yAll)] + 0.1*[-1,1]*range(yAll);

% Plot some pine trees for fun
xTree = 2.5*[0.8, 1.9, 2.8, 4.2, 8.2];
hTree = 3.5*[0.8, 1.3, 1.0, 0.7, 1.2];
for i=1:length(xTree)
    drawPineTree(xTree(i),groundModel(xTree(i)), hTree(i));
end

% Plot curves
hGnd = plot(xGround, yGround,...
    'LineWidth',groundWidth,'Color',groundColor);
hTraj = plot(xFull,yFull,...
    'LineWidth',lineWidth,'Color',lineColor);


set(gca,'XLim',xBnd);
set(gca,'YLim',yBnd);
axis equal
xlabel('horizontal position (m)')
ylabel('vertical position (m)')

% Moving parts of the animation:
frameRate = 1/20;
tFrame = tFull(1):frameRate:tFull(end);
for frame = 1:length(tFrame)
    
    % Interpolate to find the state now:
    zFrame = interp1(tFull',zFull',tFrame(frame));
    xFrame = zFrame(1);
    yFrame = zFrame(2);
    qFrame = zFrame(3);
    dqFrame = zFrame(6);
    
    % Draw the catapult
    if frame == 1
    [hArm, hTip] = drawCatapult(qFrame,dqFrame,param);
    
    % Plot the projectile:
    hProjectile = plot(xFrame,yFrame,'.','MarkerSize',40,'Color',markerColor);
    else
       
        % Position of the tip of the catapult
        zTip = getProjectileState([qFrame;dqFrame],param);
        xTip = zTip(1);
        yTip = zTip(2);
        
        % Update plot handles
        set(hArm,...
            'XData',[xTip,param.xCatapult],...
            'YData',[yTip,param.yCatapult]);
        set(hTip,'XData',xTip,'YData',yTip);
        set(hProjectile,'XData',xFrame,'YData',yFrame);
        
    end
    
    % Delay for plotting:
    drawnow;
    pause(frameRate);   % This is highly inaccurate*
    
end

