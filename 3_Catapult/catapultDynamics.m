function dz = catapultDynamics(z,param)
% dz = catapultDynamics(z,param)
%
% This function computes the dynamics of a simple torqsion catapult, in
% first-order form. An arm angle of zero corresponds to the arm pointing
% straight out in the positive vertical direction
%
% INPUTS:
%   z = [q; dq] = state
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
%   - catapult is powered by a simple torsion spring
%   - carapult arm is a slendar rod
%   - projecile is a point mass
%  

% Unpack the state
q = z(1,:);
dq = z(2,:);

% Compute the rotational inertia of the arm:
r = param.armLength;
m1 = param.projectileMass;
m2 = param.armMass;
Inertia_projectile = m1*r*r;    %Point mass
Inertia_arm = (1/3)*m2*r*r;    %slendar rod 
Inertia = Inertia_arm + Inertia_projectile;

% Compute the net torque acting on the arm:
g = param.gravity;
k = param.springConstant;
q0 = param.springRestAngle;
Torque_gravity = g*(m1*r + 0.5*m2*r)*sin(q);
Torque_spring = -k*(q-q0);
Torque = Torque_gravity + Torque_spring;

% Angular acceleration, compute by torque balance
ddq = Torque./Inertia;

% Combine into a single vector for output:
dz = [dq;ddq];

end