function [E,T,U] = energy(z,p)
%  [E,T,U] = energy(z,p)
%
% This function computes the dynamics of a disk rolling on a see-saw
%
% INPUTS:
%   z(1,:) = absolute angle of the see-saw
%   z(2,:) = absolute angle of the disk
%   z(3,:) = absolute angle rate of the see-saw
%   z(4,:) = absolute angle rate of the disk
%   p.m = mass of the disk
%   p.I1 = moment of inertia of the see-saw about CoM == pivot
%   p.I2 = moment of inertia of the disk about its CoM
%   P.R = radius of the disk
%   P.g = gravity acceleration
% 
% OUTPUTS:
%   E = total mechanical energy
%   T = kinetic energy
%   U = potential energy    
%

q1 = z(1,:);
q2 = z(2,:);
dq1 = z(3,:);
dq2 = z(4,:);

[E,T,U] = autoGen_energy(q1,q2,dq1,dq2,p.m,p.I1,p.I2,p.R,p.g);

end