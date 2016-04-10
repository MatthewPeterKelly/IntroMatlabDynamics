function zProjectile = getProjectileState(zCatapult, param)
% zProjectile = getProjectileState(zCatapult, param)
%
% Gets the state of the projectile, given the state of the catapult. This
% is used for computing the initial state for the projectile simulation.
%
% INPUTS:
%   zCatapult = [q; dq]
%   param = struct of catapult parameters
% 
% OUTPUTS:
%   zProjectile = [x;y;dx;dy]
%

x0 = param.xCatapult;
y0 = param.yCatapult;
R = param.armLength;

q = zCatapult(1,:);
dq = zCatapult(2,:);

x = x0 - R*sin(q);
y = y0 + R*cos(q);
dx = -dq.*R.*cos(q);
dy = -dq.*R.*sin(q);

zProjectile = [x;y;dx;dy];

end