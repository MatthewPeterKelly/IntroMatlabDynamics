function dz = projectileDynamics(z,param)
% dz = projectileDynamics(z,param)
%
% This function computes the dynamics of a the projectile after it has left
% the catapult and is flying through the air
%
% INPUTS:
%   z = [x; y; dx; dy] = state
%       x = projectile horizontal position
%       y = projectile vertical position
%       dx = projectile horizontal velocity
%       dy = projectile vertical velocity
%
%   param = struct of parameters
%       .armMass
%       .projectileMass
%       .armLength
%       .gravity
%       .springConstant
%       .springRestAngle
%       .quadraticAirDrag
%
% OUTPUTS:
%   dz = [dz/dt] = [dq;ddq] = derivative of state
%
% NOTES:
%   - projectile is a point mass
%   - constant gravity force
%   - quadratic air drag
%  

% Unpack the state
x = z(1,:);
y = z(2,:);
dx = z(3,:);
dy = z(4,:);

% Drag forces
speed = sqrt(dx.*dx + dy.*dy);
Cd = param.quadraticAirDrag;
Fx_drag = -Cd*speed*dx;
Fy_drag = -Cd*speed*dy;

% Gravity force:
m = param.projectileMass;
g = param.gravity;
Fy_gravity = -m*g;

% Dynamics  F=ma  -->   a = F/m
ddx = Fx_drag/m;
ddy = (Fy_drag+Fy_gravity)/m;

% Pack up derivatives:
dz = [dx;dy;ddx;ddy];

end