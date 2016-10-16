
%% Load case
%
% Airfoil rib compressed on its upper and lower edges. Stiff design.

%% Optimization type
pCon.type = 1;                                  % 1: compliance optimization
                                                % 2: compliant mechanism
%% Dimensions
pCon.Lx=4;                                      % Problem domain length
pCon.Ly=1;                                      % Problem domain height


%% Load
pCon.P = 1;


%% Geometric parameters
R = 6*pCon.Ly;
dx = 3*pCon.Lx/8;
dy = -11*pCon.Ly/2;

a = pCon.Ly/2 - R + sqrt(R^2-dx^2);
b = pCon.Ly/2 - R + sqrt(R^2-(pCon.Lx-dx)^2);

R1 = pCon.Ly/4;
R2 = pCon.Ly/4;
R3 = pCon.Ly/6;

xs1 = ((pCon.Lx/4)+R1+(((pCon.Lx/4)-R1-R2)/2)-(pCon.Lx/80));
xs2 = ((pCon.Lx/4)+R1+(((pCon.Lx/4)-R1-R2)/2)+(pCon.Lx/80));
ys1 = -2*pCon.Ly/5;
ys2 = 2*pCon.Ly/5;

%% Boundary conditions
% Lower edge compression
pCon.pLoad(1).x = [0 ; -pCon.Ly/2];
pCon.pLoad(1).F = [-pCon.P ; pCon.P];
pCon.pLoad(2).x = [dx ; -pCon.Ly/2];
pCon.pLoad(2).F = [-pCon.P ; pCon.P];
pCon.pLoad(3).x = [2*dx ; -pCon.Ly/2];
pCon.pLoad(3).F = [-pCon.P ; pCon.P];
pCon.pLoad(4).x = [pCon.Lx ; -pCon.Ly/2];
pCon.pLoad(4).F = [-pCon.P ; pCon.P];
% Upper edge compression
pCon.pLoad(5).x = [0 ; a];
pCon.pLoad(5).F = [pCon.P ; -pCon.P];
pCon.pLoad(6).x = [dx ; pCon.Ly/2];
pCon.pLoad(6).F = [pCon.P ; -pCon.P];
pCon.pLoad(7).x = [2*dx ; a];
pCon.pLoad(7).F = [pCon.P ; -pCon.P];
pCon.pLoad(8).x = [pCon.Lx ; b];
pCon.pLoad(8).F = [pCon.P ; -pCon.P];
%pCon.loads = [8 ; 0];

delta = 0.000001;
pCon.lbc(1).x = [0,-pCon.Ly/2+delta;                % Begining point of the line
                 0,a-delta];                    % End point of the line
pCon.lbc(1).type = 1;                         % Number of imposed displacements (1: clamp, 2: locking)
pCon.lbc(1).u = @(x) [0;0];                   % Imposed displacement law

% pCon.lbc(2).x = [pCon.Lx,-pCon.Ly/2;                % Begining point of the line
%                  pCon.Lx,pCon.Ly/2];                % End point of the line
% pCon.lbc(2).type = 1;                         % Number of imposed displacements (1: clamp, 2: locking)
% pCon.lbc(2).u = @(x) [0;0];                   % Imposed displacement law

%% Domain definition
pCon.domains(1).f = @(x,y) boundaryCurveAirfoil(x,y,R,dx,dy);  

%% Holes definition


% Spar
pCon.holes(1).type = 1;
pCon.holes(1).x0 = [(xs1+xs2)/2 ; (ys1+ys2)/2];
pCon.holes(1).l = [(xs2-xs1) ; (ys2-ys1)];
pCon.holes(1).coverBoundary = true;

% First circle
pCon.holes(2).type = 2;
pCon.holes(2).x0 = [pCon.Lx/4 ; -pCon.Ly/10];
pCon.holes(2).r = R1;
pCon.holes(2).coverBoundary = false;

% Second circle
pCon.holes(3).type = 2;
pCon.holes(3).x0 = [pCon.Lx/2 ; -pCon.Ly/10];
pCon.holes(3).r = R2;
pCon.holes(3).coverBoundary = false;

% Third circle
pCon.holes(4).type = 2;
pCon.holes(4).x0 = [3*pCon.Lx/4 ; -pCon.Ly/10];
pCon.holes(4).r = R3;
pCon.holes(4).coverBoundary = false;

% Fourth circle
% pCon.holes(5).type = 2;
% pCon.holes(5).x0 = [3*pCon.Lx/4 ; 0];
% pCon.holes(5).r = pCon.Ly/8;
% pCon.holes(5).coverBoundary = false;
