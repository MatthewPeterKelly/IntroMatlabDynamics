function drawCatapult(zCatapult,param)
% drawCatapult(zCatapult,param)
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
%       .quadraticAirDrag
%       .initialAngle  --  (rad) measured from pos. vert. axis
%       .launchAngle  -- (rad) measured from pos. vert. axis.
%       .xCatapult -- horizontal position of catapult axle
%       .yCatapult -- height of the catapult axle above ground

%
% OUTPUTS:
%   A pretty picture of the catapult
%


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                    figure out where stuff goes                          %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Position of the axle of the catapult
x0 = param.xCatapult;
y0 = param.yCatapult;

% Position of the tip of the arm
zProjectile = getProjectileState(zCatapult, param);
xP = zProjectile(1);
yP = zProjectile(2);
 
% Compute the sequence of positions swept by the arm:
qSweep = linspace(param.launchAngle, param.initialAngle, 25);
dqSweep = zeros(size(qSweep));
zSweep = getProjectileState([qSweep; dqSweep],param);
xSweep = zSweep(1,:);
ySweep = zSweep(2,:);

% final hard stop:
xS1 = 0.6*xSweep(1) + 0.4*x0;
yS1 = 0.6*ySweep(1) + 0.4*y0;

% initial hard stop:
xS2 = 1.2*xSweep(end) -0.2*x0;
yS2 = 1.2*ySweep(end) -0.2*y0;

% Position where the catapult support hits the ground:
xG1 = x0;
yG1 = groundModel(xG1);

xG2 = xS2;
yG2 = groundModel(xG2);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     now start drawing stuff                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Support structure:
plot([x0,xG1],[y0,yG1],'k-','LineWidth',8);
plot([xS2,xG2],[yS2,yG2],'k-','LineWidth',8);
plot([x0,xS2],[y0,yS2],'k-','LineWidth',6);
plot(x0,y0,'k.','MarkerSize',40);
plot(xS2,yS2,'k.','MarkerSize',40);

% upper hard stop:
plot([x0,xS1],[y0,yS1],'k-','LineWidth',6);
plot(xS1,yS1,'k.','MarkerSize',40);

% path swept by the arm:
plot(xSweep,ySweep,'k--','LineWidth',1);

% plot the arm itself:
plot([x0,xP],[y0,yP],'LineWidth',3,'Color',0.4*[1,1,1]);
plot(xP,yP,'ko','LineWidth',5,'MarkerSize',15,'Color',[0.2,0.2,0.8]);

end



