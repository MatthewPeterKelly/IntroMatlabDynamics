function y = groundModel(x)
% y = groundModel(x)
%
% Returns the height y of the ground at a horizontal position x
%

y = sin(0.2*x) + 0.2*cos(0.1*x).*sin(x) +0.4*cos(0.3*x);

end