% MAIN  --  Plotting
%
% This script will demonstrate how to use a variety of the basic plotting
% features of Matlab, and how to use different kinds of functions.
%


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           PART ONE                                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Create 100 points uniformly spaced between 0 and 1
t = linspace(0,1,100);
x = sin(4*t); 
v = 4*cos(4*t);

% Create a simple plot
figure(1); clf; 
plot(t,x);
xlabel('time (s)')
ylabel('position (m)')
title('Position vs Time')


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           PART TWO                                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% How to show two different quantities?
figure(2); clf;

% Creates a rectangular grid of plots, 
nRows = 2;
nCols = 1;

idxSelect = 1;
subplot(nRows,nCols,idxSelect);  
plot(t,x);
xlabel('time (s)')
ylabel('position (m)')
title('Position vs Time')

idxSelect = 2;
subplot(nRows,nCols,idxSelect);  
plot(t,v);
xlabel('time (s)')
ylabel('velocity (m/s)')
title('Velocity vs Time')


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           PART THREE                                    %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% What if there were two masses now, instead of one?

% Create data:
t = linspace(0,2,100);
x1 = sin(4*t); 
v1 = 4*cos(4*t);
x2 = 0.6*sin(9*t); 
v2 = 0.6*9*cos(9*t);

% Open and clear a new figure:
figure(3); clf;

% Creates a rectangular grid of plots, 
nRows = 2;
nCols = 1;

idxSelect = 1;
subplot(nRows,nCols,idxSelect); 
hold on;   %Do this to plot two lines together
plot(t,x1);
plot(t,x2);
xlabel('time (s)')
ylabel('position (m)')
title('Position vs Time')
legend('x1','x2')   %This makes a legend to see each different curve

idxSelect = 2;
subplot(nRows,nCols,idxSelect); 
hold on;   %Do this to plot two lines together
plot(t,x1);
plot(t,x2);
xlabel('time (s)')
ylabel('velocity (m/s)')
title('Velocity vs Time')
legend('v1','v2')   %This makes a legend to see each different curve



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           PART Four                                     %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Lets suppose that we want to use a single function call for both position
% and velocity:
myFun = @(t)( [sin(4*t); 4*cos(4*t)] );

z = myFun(t);   %Call an anonymous function to create a matrix z
x = z(1,:);
y = z(2,:);

% Let's plot the particle in phase space (velocity vs position)
figure(4); clf;
plot(x,y);


