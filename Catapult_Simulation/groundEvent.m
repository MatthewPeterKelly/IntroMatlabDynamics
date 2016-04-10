function [value, isTerminal, direction] = groundEvent(z)
% [value, isTerminal, direction] = groundEvent(z)
%
% Triggers an event that tells ode45 to stop the simulation when the
% projectile reaches the ground.
%
% INPUTS:
%   z = [x;y;dx;dy] = state
%       x = horizontal position of projectile
%       y = vertical position of the projectile
%
% OUTPUTS:
%   value = ode45 should trigger and event when value --> 0
%   direction = trigger only when the value function is
%       "1"  -->  trigger only when value is increasing
%       "0"  -->  trigger whenever value is zero
%       "-1" -->  trigger only when value is decreasing
%   isTerming = should ode45 stop simulation when event triggers?
%

x = z(1,:);  %horizontal position of the projecile
yP = z(2,:);   %vertical position of the projectile
yG = groundModel(x);
value = yG - yP;  %Trigger event when this is zero
direction = zeros(size(value));  %Always trigger
isTerminal = true(size(value));  %stop on event

end