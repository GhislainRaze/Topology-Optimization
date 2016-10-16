
%% Load case
%
% Airfoil rib compressed on its upper and lower edges. Compliant design
% with the edges compression as input and the spar compression as output.

%% Optimization type
pCon.type = 2;                                  % 1: compliance optimization
                                                % 2: compliant mechanism
%% Dimensions
pCon.Lx=4;                                      % Problem domain length
pCon.Ly=1;                                      % Problem domain height


%% Load
pCon.P = 1;


%% Boundary conditions
% Lower edge compression
pCon.lLoad(1).x = [0 ,-pCon.Ly/2;
                   pCon.Lx, -pCon.Ly/2];
pCon.lLoad(1).F =  @(x) [0 ; pCon.P];
% Upper edge compression
pCon.lLoad(2).x = [0 ,pCon.Ly/4;
                    pCon.Lx, pCon.Ly/4];
pCon.lLoad(2).F = @(x) [0 ; -pCon.P];
pCon.lLoad(2).param = @(x) [x;boundaryCurveAirfoil(x,0)];
pCon.lLoad(2).length = pCon.Lx;

    
pCon.loads = [2 6 ; 0 0];

pCon.lbc(1).x = [0,-pCon.Ly/2;                % Begining point of the line
                 0,pCon.Ly/2];                % End point of the line
pCon.lbc(1).type = 1;                         % Number of imposed displacements (1: clamp, 2: locking)
pCon.lbc(1).u = @(x) [0;0];                   % Imposed displacement law

pCon.lbc(2).x = [pCon.Lx,-pCon.Ly/2;                % Begining point of the line
                 pCon.Lx,pCon.Ly/2];                % End point of the line
pCon.lbc(2).type = 1;                         % Number of imposed displacements (1: clamp, 2: locking)
pCon.lbc(2).u = @(x) [0;0];                   % Imposed displacement law

%% Domain definition
pCon.domains(1).f = @(x,y) boundaryCurveAirfoil(x,y);  

%% Holes definition
% First circle
pCon.holes(1).type = 2;
pCon.holes(1).x0 = [pCon.Lx/8 ; -pCon.Ly/8];
pCon.holes(1).r = pCon.Ly/6;
pCon.holes(1).coverBoundary = false;

% Spar
pCon.holes(2).type = 1;
pCon.holes(2).x0 = [pCon.Lx/4 ; 0];
pCon.holes(2).l = [pCon.Lx/20 ; 7*pCon.Ly/8];
pCon.holes(2).coverBoundary = true;

% Second circle
pCon.holes(3).type = 2;
pCon.holes(3).x0 = [3*pCon.Lx/8 ; -pCon.Ly/8];
pCon.holes(3).r = pCon.Ly/6;
pCon.holes(3).coverBoundary = false;

% Third circle
pCon.holes(4).type = 2;
pCon.holes(4).x0 = [3*pCon.Lx/4 ; 0];
pCon.holes(4).r = pCon.Ly/8;
pCon.holes(4).coverBoundary = false;


% Spar compression load
pCon.lLoad(3).x = [pCon.holes(2).x0(1)-pCon.holes(2).l(1)/2, pCon.holes(2).x0(2)-pCon.holes(2).l(2)/2 ;
                    pCon.holes(2).x0(1)+pCon.holes(2).l(1)/2, pCon.holes(2).x0(2)-pCon.holes(2).l(2)/2 ];
pCon.lLoad(3).F =  @(x) [0 ; pCon.P];      
pCon.lLoad(4).x = [pCon.holes(2).x0(1)+pCon.holes(2).l(1)/2, pCon.holes(2).x0(2)-pCon.holes(2).l(2)/2 ;
                    pCon.holes(2).x0(1)+pCon.holes(2).l(1)/2, pCon.holes(2).x0(2)+pCon.holes(2).l(2)/2 ];
pCon.lLoad(4).F =  @(x) [-pCon.P ; 0];  
pCon.lLoad(5).x = [pCon.holes(2).x0(1)+pCon.holes(2).l(1)/2, pCon.holes(2).x0(2)+pCon.holes(2).l(2)/2 ;
                    pCon.holes(2).x0(1)-pCon.holes(2).l(1)/2, pCon.holes(2).x0(2)+pCon.holes(2).l(2)/2 ];
pCon.lLoad(5).F =  @(x) [0 ; -pCon.P];  
pCon.lLoad(6).x = [pCon.holes(2).x0(1)-pCon.holes(2).l(1)/2, pCon.holes(2).x0(2)+pCon.holes(2).l(2)/2 ;
                    pCon.holes(2).x0(1)-pCon.holes(2).l(1)/2, pCon.holes(2).x0(2)-pCon.holes(2).l(2)/2 ];
pCon.lLoad(6).F =  @(x) [pCon.P ; 0];            