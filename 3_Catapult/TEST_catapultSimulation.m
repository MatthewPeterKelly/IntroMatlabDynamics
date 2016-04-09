% TEST  --  Catapult Simulation
%
% Simulated a simple torsion catapult launching a projectile. 
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
z0 = [param.initialAngle; 0];
tSpan = [0,4];  
dynFun = @(t,z)( catapultDynamics(z,param) );
odeOpt = odeset(...
    'Event',@(t,z)( launchEvent(z,param) ),...
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
qGrid = sol.y(1,:);
dqGrid = sol.y(2,:);

%%%% Interpolat solution between grid-points
t = linspace(tGrid(1), tGrid(end), 100);  %time grid for plotting results
z = deval(sol,t);  % state at each grid point in time
q = z(1,:);  % angle of the catapult arm
dq = z(2,:); % angular rate of the catapult arm

%%%% set up for plotting
figure(1); clf; 
lineWidth = 2;
lineColor = [0.2, 0.2, 0.8];  % Blue
markerSize = 7;
markerType = 'o';
markerLine = 2;
markerColor = [0.8,0.2,0.2];  % Red

%%%% arm angle plot
subplot(2,1,1); hold on;
plot(t([1,end]),param.launchAngle*[1,1],'k--','LineWidth',1);
plot(t,q,'LineWidth',lineWidth,'Color',lineColor); 
plot(tGrid,qGrid,markerType,...
    'MarkerSize',markerSize, 'LineWidth',markerLine,'Color', markerColor);
xlabel('time (s)')
ylabel('angle (rad)')
legend('launch angle','arm angle','ode45 grid')
title('Catapult Launch Simulation')

%%%% arm angular rate plot
subplot(2,1,2); hold on;
plot(t,dq,'LineWidth',lineWidth,'Color',lineColor);  
plot(tGrid,dqGrid,markerType,...
    'MarkerSize',markerSize, 'LineWidth',markerLine,'Color', markerColor);
xlabel('time (s)')
ylabel('rate (rad/s)')
legend('arm rate','ode45 grid')
