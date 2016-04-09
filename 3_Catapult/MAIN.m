% MAIN  --  Catapult Simulation
%
% Simulated a simple torsion catapult launching a projectile. 
%
% Uses ode45 for simulation and event-detection. 
%
% Phase One:  simulation of the catapult launching the projectil
% Phase Two:  simulation of the projectile flying through the air
%


%%%% Problem parameters
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
%                  ode45 simulation of launch                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
z0Catapult = [param.initialAngle; 0];
tSpan = [0,4];  
dynCatapult = @(t,z)( catapultDynamics(z,param) );
odeOpt = odeset(...
    'Event',@(t,z)( launchEvent(z,param) ),...
    'AbsTol',1e-8,...
    'RelTol',1e-8);
solCatapult = ode45(dynCatapult,tSpan,z0Catapult,odeOpt);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                    analysis of the simulation                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%% Check to make sure that the catapult launched
zEvent = solCatapult.ye;   % state of the system at the instant of launch
if isempty(zEvent)
    error('No launch event detected. Invalid parameter set.');    
end









%%%% TODO:  Finish combined simulation and animation %%%%








%%%% Extract the grid points that ode45 actually used:
tGrid = solCatapult.x;
qGrid = solCatapult.y(1,:);
dqGrid = solCatapult.y(2,:);

%%%% Interpolat solution between grid-points
t = linspace(tGrid(1), tGrid(end), 100);  %time grid for plotting results
z = deval(solCatapult,t);  % state at each grid point in time


