% TEST  --  Projectile Simulation
%
% Simulate a projectile flying through the air with quadratic drag
%
% Uses ode45 for simulation and event-detection. 
%

%%%% parameters for the model and simulation:
param.armMass = 8;   %(kg)  arm is a slendar rod
param.projectileMass = 1;  %(kg)
param.armLength = 1.5;  %(m) 
param.gravity = 9.81;  %(m/s^2)
param.springConstant = 1000;  %(N/rad)
param.springRestAngle = 0*(pi/180);  % (rad)  measured from pos. vert. axis.
param.initialAngle = (90+30)*(pi/180);  % (rad)  measured from pos. vert. axis.
param.quadraticAirDrag = 0.1;  %(N-s^2/m^2)
param.launchAngle = 45*(pi/180);   %(rad) measured from pos. vert. axis.


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        ode45 simulation                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
z0 = [0;3;6;5];  % [x;y;dx;dy]  initial projectile state
tSpan = [0,4];  
dynFun = @(t,z)( projectileDynamics(z,param) );
odeOpt = odeset(...
    'Event',@(t,z)( groundEvent(z) ),...
    'AbsTol',1e-8,...
    'RelTol',1e-8);
sol = ode45(dynFun,tSpan,z0,odeOpt);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                    analysis of the simulation                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%% Check to make sure that the catapult launched
zEvent = sol.ye;   % state of the system at the instant of launch
if isempty(zEvent)
    error('No launch event detected. Invalid parameter set.');    
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
groundWidth = 3;
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

%%%% Draw a picture of the whole thing on a new figure:
figure(3); clf; hold on;

xBnd = [min(x), max(x)] + 0.1*[-1,1]*range(x);
xGround = linspace(xBnd(1), xBnd(2), 150);
yGround = groundModel(xGround);
yAll = [y,yGround];  %Collect all y points, used for scaling only
yBnd =  [min(yAll), max(yAll)] + 0.1*[-1,1]*range(yAll);


plot(xGround, yGround, 'LineWidth',groundWidth,'Color',groundColor); 
plot(x,y,'LineWidth',lineWidth,'Color',lineColor);  
plot(xGrid,yGrid,markerType,...
    'MarkerSize',markerSize, 'LineWidth',markerLine,'Color', markerColor);
set(gca,'XLim',xBnd);
set(gca,'YLim',yBnd);
legend('ground','trajectory','ode45 grid');
axis equal
xlabel('horizontal position (m)')
ylabel('vertical position (m)')

