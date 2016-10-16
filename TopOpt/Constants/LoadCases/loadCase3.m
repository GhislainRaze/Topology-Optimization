
%% Load case 2
%
% L-shape problem


%% Optimization type
pCon.type = 1;                                  % 1: compliance optimization
                                                % 2: compliant mechanism
%% Dimensions
pCon.Lx=2;                                      % Problem domain length
pCon.Ly=2;                                      % Problem domain height


%% Load
pCon.P = 1;

%% Additionnal domains
shift = 1e-6;
pCon.holes(1).type = 1;
pCon.holes(1).l = [pCon.Lx/2 ; pCon.Ly/2];
pCon.holes(1).x0 = [3*pCon.Lx/4+shift;pCon.Ly/4+shift];
pCon.holes(1).coverBoundary = false;


%% Boundary conditions
pCon.pLoad(1).x = [pCon.Lx ;0];
pCon.pLoad(1).F = [0 ; -pCon.P];

pCon.lbc(1).x = [0,pCon.Ly/2;                % Begining point of the line
                 pCon.Lx/2,pCon.Ly/2];                % End point of the line
pCon.lbc(1).type = 1;                         % Number of imposed displacements (1: clamp, 2: locking)
pCon.lbc(1).u = @(x) [0;0];                   % Imposed displacement law
    