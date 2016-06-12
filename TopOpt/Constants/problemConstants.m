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
% * _pCon.vol_: the total volume of the structure (regardless of the mass
% distribution)
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
    pCon.holes = [];

    
    %% Natural boundary conditions
    
    % Body forces
    pCon.b=[0;0]; 
    
    % Line forces
    pCon.lLoad=[];
    
    % Point forces
    pCon.pLoad = [];
    
    %% Essential boundary conditions
    
    % Line boundary conditions
    pCon.lbc = [];
    
    
    % Point boundary conditions
    pCon.pbc = [];

    %% Load case determination
    eval(loadCase);
    
    %% Volume computation
    pCon.vol = pCon.Lx*pCon.Ly;
    for i = 1 : length(pCon.holes)
        if pCon.holes(i).type == 1
            pCon.vol = pCon.vol - pCon.holes(i).l(1)*pCon.holes(i).l(2);
        elseif pCon.holes(i).type == 2
            pCon.vol = pCon.vol - pi*pCon.holes(i).r^2;
        end
    end
    
    %% Boundary conditions additional parameters
    pCon = boundaryConditionsParameters(pCon);

end