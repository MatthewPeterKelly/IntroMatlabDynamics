% MAIN  --  Catapult Simulation
%
% Simulated a simple torsion catapult launching a projectile. 
%
% Uses ode45 for simulation and event-detection. 
%
% Phase One:  simulation of the catapult launching the projectil
% Phase Two:  simulation of the projectile flying through the air
%



%%%% parameters for the catapult and projectile model:
param.armMass = 8;   %(kg)  arm is a slendar rod
param.projectileMass = 1;  %(kg)
param.armLength = 1.5;  %(m) 
param.gravity = 9.81;  %(m/s^2)
param.springConstant = 1000;  %(N/rad)
param.springRestAngle = 0*(pi/180);  % (rad)  measured from pos. vert. axis.
param.initialAngle = (90+30)*(pi/180);  % (rad)  measured from pos. vert. axis.
param.quadraticAirDrag = 0.1;  %(N-s^2/m^2)

%%%% parameters for the simulation
param.launchAngle = 45*(pi/180);   %(rad) measured from pos. vert. axis.


%%%% set up catapult simulation
z0 = [param.initialAngle; 0];
tSpan = [0,1];
dynFun = @(t,z)( catapultDynamics(z,param) );
odeOpt = odeset(...
    'Event',@(t,z)( launchEvent(z,param) ),...
    'AbsTol',1e-8,...
    'RelTol',1e-8);
sol = ode45(dynFun,tSpan,z0,odeOpt);





