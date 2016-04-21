%%%% Derive Equations of Motion
%
% "See-Saw" Disk problem
%
% This system consists of two components:
%
% 1) The see-saw, which is modeled as a slendar rod that rotates about its
% center of mass.
%
% 2) A uniform disk, which rolls without slip along the see-saw
%
% m = mass of disk
% R = radius of the disk
% 
% I1 = moment of inertia of the rod, about its CoM
% I2 = moment of inertia of the disk, about its CoM
%
% q1 = absolute angle of the rod, measured from horizontal
% q2 = absolute angle of the disk, measured from horizontal
%
%

clc; clear;

syms m I1 I2 R g 'real'
syms q1 q2 dq1 dq2 ddq1 ddq2 'real'

%%%% State vector and derivative:
z = [q1;q2;dq1;dq2];
dz = [dq1;dq2;ddq1;ddq2];
derivative = @(in)( jacobian(in,z)*dz );  %Chain rule!

%%%% Unit vectors:

% Inertial reference frame
e1 = sym([1;0;0]);   % horizontal
e2 = sym([0;1;0]);   % vertical
e3 = sym([0;0;1]);   % out of page

% Frame fixed to the see-saw
a1 = cos(q1)*e1 + sin(q1)*e2;  %pointing along rod
a3 = e3;  %out of the page
a2 = cross(a3,a1);  %pointing normal to to the rod;


%%%% Rolling without slip:
x = -(q2-q1)*R;   %distance of contact point along rod

%%%% Kinematics 
% o = origin, pivot point for rod (see-saw)
% c = contact point between disk and rod
% p = center of the disk

% Position vectors:
Rco = x*a1;   
Rpc = R*a2;
Rpo = Rpc + Rco;

% Derivatives, in the inertial reference frame
dRpo = derivative(Rpo);
ddRpo = derivative(dRpo);

%%%% Forces:
% The only force that exerts an external moment about the pivot of the rod
% is the weight of the disk. 

Fg = -m*g*e2;   %Weight of the disk


%%%% Angular Momentum Balance:
% Entire System, about pivot point (o)
sumTorques1 = cross(Rpo,Fg);
sumInertial1 = cross(Rpo,m*ddRpo) + ddq1*I1*e3 + ddq2*I2*e3;
eqn1 = simplify(dot(sumTorques1-sumInertial1,e3));


%%%% Angular Momentum Balance:
% Disk Only, about the contact point (c)
sumTorques2 = cross(Rpc,Fg);
sumInertial2 = cross(Rpc,m*ddRpo) + ddq2*I2*e3;
eqn2 = simplify(dot(sumTorques2-sumInertial2,e3));


%%%% Solve Equations of motion:
eqns = [eqn1;eqn2];
vars = [ddq1;ddq2];
[MM,ff] = equationsToMatrix(eqns,vars);
soln = simplify(MM\ff);


%%%% Display the solution:
fprintf('\n======================= \nddq1 = \n')
pretty(soln(1));
fprintf('\n======================= \nddq2 = \n')
pretty(soln(2));


%%%% Check energy conservation:
U = m*g*dot(Rpo,e2);
T = (1/2)*I1*dq1^2 + (1/2)*I2*dq2^2 + (1/2)*m*dot(dRpo,dRpo);
E = simplify(U + T);


%%%% Write out function files:
matlabFunction(soln(1), soln(2),...
    'file','autoGen_dynamics.m',...
    'vars',{q1,q2,dq1,dq2,m,I1,I2,R,g},...
    'outputs',{'ddq1','ddq2'});
matlabFunction(E,T,U,...
    'file','autoGen_energy.m',...
    'vars',{q1,q2,dq1,dq2,m,I1,I2,R,g},...
    'outputs',{'E','T','U'});

