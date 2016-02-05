% MAIN  --  Script
%
% A simple introduction to scripts in Matlab.

%% What is a script?
%
% A script is a matlab file (filename.m) that runs as if it was directly
% entered into the command prompt. As a result, all variables that you
% create remain in the workspace after the script is done running.
%
% You can think of a script as a simple tool to write down a list of
% commands that you want to run one after the other

% This is a command for writing text to the command prompt:
disp('Hello World!');  

%% How to run a script?
%
% You can run a script by typing the file name in the command prompt:
% >> MAIN_0_script
% This assumes that you have the script file on the Matlab path.
%
% You can click on the green play button that says "Run" on the top of the
% window, as long as the script is selected in the editor.
%
% On Windows (and Linux) you can press F5 to run a script as well. There is
% a similar command on Mac, but I don't have one to test it on now.

%%  A script can do basic operations and print them to the command prompt

% The result will be assigned to the variable "ans"
2+3

% We can now use "ans".
value = ans + 4

% Here were are going to modify value, but not print the result. 
value = value + 1;   %Notice the semicolon!

%% A script can call a function:

% Create 100 points, evenly spaced between zero and five
t = linspace(0,5,100);

% Call a function. Suppress output.
x = sin(t);

%% We can also plot using a script:
figure(1); clf;  % Open Figure 1 and clear it
plot(t,x);  %Plot x vs t





