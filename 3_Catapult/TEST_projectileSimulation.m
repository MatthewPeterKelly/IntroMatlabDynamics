% TEST  --  Projectile Simulation
%
% Simulate a projectile flying through the air with quadratic drag
%
% Uses ode45 for simulation and event-detection. 
%

%%%% parameters for the model and simulation:
param.armMass = 10;   %(kg)  arm is a slendar rod
param.projectileMass = 2;  %(kg)
param.armLength = 2.5;  %(m) 
param.gravity = 9.81;  %(m/s^2)
param.springConstant = 3000;  %(N/rad)
param.springRestAngle = 0*(pi/180);  % (rad)  measured from pos. vert. axis.
param.initialAngle = (90+30)*(pi/180);  % (rad)  measured from pos. vert. axis.
param.quadraticAirDrag = 0.1;  %(N-s^2/m^2)
param.launchAngle = 45*(pi/180);   %(rad) measured from pos. vert. axis.
param.xCatapult = 0;  %(m)  horizontal position of catapult axle
param.yCatapult = 2;  %(m)  height of the catapult axle above ground

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        ode45 simulation                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
q0 = param.launchAngle;
dq0 = -4;  %Angular rate of catapult during launch
z0Catapult = [q0;dq0];
z0Projectile = getProjectileState(z0Catapult, param);

tSpan = [0,4];  
dynFun = @(t,z)( projectileDynamics(z,param) );
odeOpt = odeset(...
    'Event',@(t,z)( groundEvent(z) ),...
    'AbsTol',1e-8,...
    'RelTol',1e-8);
sol = ode45(dynFun,tSpan,z0Projectile,odeOpt);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                    analysis of the simulation                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%% Check to make sure that the catapult launched
zEvent = sol.ye;   % state of the system at the instant of launch
if isempty(zEvent)
    error('No ground event detected. Invalid parameter set.');    
end

%%%% Extract the grid points that ode45 actually used:
tGrid = sol.x;
xGrid = sol.y(1,:);
yGrid = sol.y(2,:);
dxGrid = sol.y(3,:);
dyGrid = sol.y(4,:);

%%%% Interpolat solution between grid-points
t = linspace(tGrid(1), tGrid(end), 100);  %time grid for plotting results
z = deval(sol,t);  % state at each grid point in time
x = z(1,:);  
y = z(2,:);
dx = z(3,:);
dy = z(4,:);
yG = groundModel(x);

%%%% set up for plotting
figure(2); clf; 
lineWidth = 2;
lineColor = [0.2, 0.2, 0.8];  % Blue
markerSize = 7;
markerType = 'o';
markerLine = 2;
markerColor = [0.8,0.2,0.2];  % Red
groundWidth = 5;
groundColor = [77,38,0]/255;

%%%% horizontal position
subplot(2,2,1); hold on;
plot(t,x,'LineWidth',lineWidth,'Color',lineColor); 
plot(tGrid,xGrid,markerType,...
    'MarkerSize',markerSize, 'LineWidth',markerLine,'Color', markerColor);
xlabel('time (s)')
ylabel('horizontal position (m)')
legend('arm angle','ode45 grid')

%%%% horizontal velocity
subplot(2,2,3); hold on;
plot(t,dx,'LineWidth',lineWidth,'Color',lineColor);  
plot(tGrid,dxGrid,markerType,...
    'MarkerSize',markerSize, 'LineWidth',markerLine,'Color', markerColor);
xlabel('time (s)')
ylabel('rate (rad/s)')
legend('arm rate','ode45 grid')


%%%% vertical position
subplot(2,2,2); hold on;
plot(t,yG,'LineWidth',groundWidth,'Color',groundColor); 
plot(t,y,'LineWidth',lineWidth,'Color',lineColor); 
plot(tGrid,yGrid,markerType,...
    'MarkerSize',markerSize, 'LineWidth',markerLine,'Color', markerColor);
xlabel('time (s)')
ylabel('vertical position (m)')
legend('ground','arm angle','ode45 grid')

%%%% vertical velocity
subplot(2,2,4); hold on;
plot(t,dy,'LineWidth',lineWidth,'Color',lineColor);  
plot(tGrid,dyGrid,markerType,...
    'MarkerSize',markerSize, 'LineWidth',markerLine,'Color', markerColor);
xlabel('time (s)')
ylabel('rate (rad/s)')
legend('arm rate','ode45 grid')



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         draw a picuture!                                %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
figure(3); clf; hold on;

% Get the scaling for the graphics
xBnd = [min(x), max(x)] + 0.2*[-1,1]*range(x);
xGround = linspace(xBnd(1), xBnd(2), 150);
yGround = groundModel(xGround);
yAll = [y,yGround];  %Collect all y points, used for scaling only
yBnd =  [min(yAll), max(yAll)] + 0.1*[-1,1]*range(yAll);

% Plot some pine trees for fun
xTree = [0.8, 1.9, 2.8, 4.2, 8.2];
hTree = 1.5*[0.8, 1.3, 1.0, 0.7, 1.2];
for i=1:length(xTree)
    drawPineTree(xTree(i),groundModel(xTree(i)), hTree(i));
end

% Draw the catapult
drawCatapult(z0Catapult(1),z0Catapult(2),param);

% Plot curves
hGnd = plot(xGround, yGround,...
    'LineWidth',groundWidth,'Color',groundColor); 
hTraj = plot(x,y,...
    'LineWidth',lineWidth,'Color',lineColor);  
hGrid = plot(xGrid,yGrid,markerType,...
    'MarkerSize',markerSize, 'LineWidth',markerLine,'Color', markerColor);

% Scale and label the axis
set(gca,'XLim',xBnd);
set(gca,'YLim',yBnd);
legend([hGnd;hTraj;hGrid],{'ground';'trajectory';'ode45 grid'});
axis equal
xlabel('horizontal position (m)')
ylabel('vertical position (m)')

