
%% Load case 1
%
% Cantilever beam with dimensions 1x2 and a unit downwards point load at
% (1,0)

%% Dimensions
pCon.Lx=1;                                      % Problem domain length
pCon.Ly=2;                                      % Problem domain height


%% Load
pCon.P = 1;


%% Boundary conditions
pCon.pLoad(1).x = [pCon.Lx ;0];
pCon.pLoad(1).F = [0 ; -pCon.P];

pCon.lbc(1).x = [0,-pCon.Ly/2;                % Begining point of the line
                 0,pCon.Ly/2];                % End point of the line
pCon.lbc(1).type = 1;                         % Number of imposed displacements (1: clamp, 2: locking)
pCon.lbc(1).u = @(x) [0;0];                   % Imposed displacement law