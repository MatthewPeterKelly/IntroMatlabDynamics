function dz = projectileDynamics(z,param)
% dz = projectileDynamics(z,param)
%
% This function computes the dynamics of a simple torqsion catapult, in
% first-order form.
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
I_projectile = m1*r*r;
I_arm = (1/3)*m2*r*r;
Inertia = I_arm + I_projectile;

% Compute the net torque acting on the arm:
g = param.gravity;
k = param.springConstant;
q0 = param.springRestAngle;
T_gravity = -g*(m1*r + 0.5*m2*r)*sin(q);
T_spring = -k*(q-q0);
Torque = T_gravity + T_spring;

% Angular acceleration, compute by torque balance
ddq = Torque./Inertia;

% Combine into a single vector for output:
dz = [dq;ddq];

end