%% Problem Constants
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% This function creates the set of problem constants _pCon_. The problem
% constants are:
%
% * _pCon.Lx_: the problem domain length
% * _pCon.Ly_: the problem domain height
% * _pCon.E_: the material Young modulus
% * _pCon.nu_: the material Poisson ratio
% * _pCon.D_: the 2D Hooke tensor
%
% This function also defines the boundary conditions.
function pCon = problemConstants()

    %Problem definition constants
    pCon.Lx=1;                                      % Problem domain length
    pCon.Ly=2;                                      % Problem domain height
    pCon.E=10;                                       % Elasticiy modulus
    pCon.nu=0.3;                                    % Poisson's ratio
    pCon.D=pCon.E/(1-pCon.nu^2)*[1       pCon.nu 0; % Plain stress elasticy matrix
                                 pCon.nu 1       0; 
                                 0       0       (1-pCon.nu)/2];
    pCon.P = 0.1;
    pCon.I = pCon.Ly^3/12;
    
    %% Natural boundary conditions
    
    % Body forces
    pCon.b=[0;0]; 
    
    % Line forces
    pCon.lLoad=[];
    pCon.lLoad(1).x = [pCon.Lx, -pCon.Ly/10;        % Begining point of the line
                       pCon.Lx, pCon.Ly/10];        % End point of the line
    pCon.lLoad(1).F = @(x) [0;-pCon.P];      % Associated force law
    
    % Point forces
    pCon.pLoad = [];
%     pCon.pLoad(1).x = [pCon.Lx ; -pCon.Ly/2];
%     pCon.pLoad(1).F = [0 ; -pCon.P];
    %% Essential boundary conditions
    
    % Line boundary conditions
    pCon.lbc = [];
    pCon.lbc(1).x = [0,-pCon.Ly/2;                % Begining point of the line
                     0,pCon.Ly/2];                % End point of the line
    pCon.lbc(1).type = 1;                         % Number of imposed displacements (1: clamp, 2: locking)
    pCon.lbc(1).u = @(x) [0;0];                   % Imposed displacement law
    
    % Point boundary conditions
    pCon.pbc = [];

    %% Boundary conditions additional parameters
    pCon = boundaryConditionsParameters(pCon);

end