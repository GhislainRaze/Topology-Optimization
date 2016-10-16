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
% constants are
%
% * _pCon.Lx_: the problem domain length
% * _pCon.Ly_: the problem domain height
% * _pCon.E_: the material Young modulus
% * _pCon.nu_: the material Poisson ratio
% * _pCon.D_: the 2D Hooke tensor
%
% * _pCon.holes_: the holes in the structure. There is no node inside a
% hole.
% * _pCon.domains_: level-set functions describing the designable domain
% * _pCon.vol_: the total volume of the structure (regardless of the mass
% distribution). It is computed once the discretization has been built.
%
% This function also defines the boundary conditions
%
% * _pCon.b_: body forces
% * _pCon.lLoad_: line loads
% * _pCon.pLoad_: point loads
% * _pCon.lbc_: line essential boundary conditions
% * _pCon.pbc_: point essential boundary conditions

function problemConstants(loadCase)
    
    global pCon
    
    %Problem definition constants
    pCon.E=1;                                       % Elasticiy modulus
    pCon.nu=0.3;                                    % Poisson's ratio
    pCon.D=1/(1-pCon.nu^2)*[1       pCon.nu 0; % Plane stress stiffness matrix
                                 pCon.nu 1       0; 
                                 0       0       (1-pCon.nu)/2];
    
    %
    pCon.domains = [];
                             
    %
    pCon.holes = [];

    %
    pCon.filledRegions = [];
    %% Natural boundary conditions
    
    % Body forces
    pCon.b=[0;0]; 
    
    % Line forces
    pCon.lLoad=[];
    
    % Point forces
    pCon.pLoad = [];
    
    % Load groups description
    pCon.loads = [];
    
    %% Essential boundary conditions
    
    % Line boundary conditions
    pCon.lbc = [];
    
    
    % Point boundary conditions
    pCon.pbc = [];

    %% Load case determination
    eval(loadCase);
    
    %% Volume computation
    pCon.vol = 0;
    
    %% Boundary conditions additional parameters
    pCon = boundaryConditionsParameters(pCon);

end