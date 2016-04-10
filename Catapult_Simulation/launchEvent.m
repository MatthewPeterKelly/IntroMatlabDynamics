function [value, isTerminal, direction] = launchEvent(z,param)
% [value, isTerminal, direction] = launchEvent(z,param)
%
% Triggers an event that tells ode45 to stop the simulation when the
% catapult reaches the launch angle.
%
% INPUTS:
%   z = [q; dq] = state
%       q = arm angle, (rad) measured from pos. horiz. axis
%   param = struct of parameters
%       .armMass
%       .projectileMass
%       .armLength
%       .gravity
%       .springConstant
%       .springRestAngle
%       .quadraticAirDrag
%       .launchAngle      %(rad) measured from pos. horiz. axis.
%
% OUTPUTS:
%   value = ode45 should trigger and event when value --> 0
%   direction = trigger only when the value function is
%       "1"  -->  trigger only when value is increasing
%       "0"  -->  trigger whenever value is zero
%       "-1" -->  trigger only when value is decreasing
%   isTerming = should ode45 stop simulation when event triggers?
%

qArm = z(1,:);   % Angle of the catapult
qLaunch = param.launchAngle;  % Critical angle to fire at
value = qArm - qLaunch;  %Trigger event when this is zero
direction = zeros(size(value));  %Always trigger
isTerminal = true(size(value));  %stop on event

end