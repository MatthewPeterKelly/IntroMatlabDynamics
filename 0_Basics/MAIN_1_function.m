function output = MAIN_1_function(input)
% output = MAIN_1_function(input)
%
% This is an example that describes what a function is, and how to use it.
%
% INPUTS:
%   input = scalar = some number to manipulate
% 
% OUTPUTS:
%   output = scalar = another number
%
% NOTES:
% 
%   The name of a function should always match its file name. 
%
%   There can be many inputs and many outputs. In this example, we only
%   have one of each.
%
%   A function is different than a script in many ways. The most important
%   is that any variables that are created inside a function are removed
%   from the workspace after the function finishes.
%
%   The comments that are immediately after the function, including this
%   sentence, form the help file for your function. Try it at the command
%   prompt:
%   >> help MAIN_1_function
%
%   I use this help file format for all of my functions. The first line is
%   an exact copy of the function declaration, then the next part is a
%   sentence that describes what the function does. Then I list the INPUTS,
%   the OUTPUTS, and any NOTES about the function.
%
%   A function can be called from inside another function, a script, or a
%   command prompt. If it is called by pressing F5 or the "Run" button,
%   then Matlab will attempt to run it without any arguments.

% A function can print things to command prompt as well:
disp('Hello World!');

% Check to make sure that we've received enough inputs:
if nargin < 1
    error('Not enough inputs!')
end

% Do basic operations:
output = input + 5;

% We can check the number of inputs that were passed. This line is
% complicated. The [ ] brackets are forming an array of strings. The first
% string is explicitly written out. The second string is created by a call
% to a matlab built-in function called num2str that converts a number into
% a string. The argument to that function is another matlab function called
% nargin, which returns the number of arguments passed to this function.
disp(['Number of inputs: ' num2str(nargin)]);


% We can call functions inside of functions.
timeMs = getTimeMs();
disp(['Time since Matlab started: ' num2str(timeMs) ' ms']);

% We can also define a function inside of another function:
a = 2.3;
b = 5.4;
C = 6.5;
D = 0.5;
    function val = addVals(a,b)
        % Notice that there are two ways to pass variables to this type of
        % function. a and b are passed as arguments to the function. Since
        % C and D are not defined inside of this function, Matlab knows to
        % look for them in this functions parent. 
        val = a + b + C + D;        
    end

% There are even more kinds of functions. An anonymous function is one that
% does not have a name. It is defined using the @ symbol. The resulting
% variable is called a function handle.
myFunHandle = @(t)( t.^2 );

% We can call a function handle just like a regular function
myFunHandle(3.2)

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function timeMs = getTimeMs()
%
% This is a sub-function. It can only be called from inside of this
% function. In other words, it cannot be called from the command prompt.

timeMs = 1000*cputime;

end


