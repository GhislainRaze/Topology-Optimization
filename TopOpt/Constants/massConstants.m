%% Mass Constants
%
%
% This functions takes a set of problem constants _pCon_ and a set of
% discretization constants _mCon_ as arguments. It is used to define a
% series of mass distribution parameters :
%
% * _mmCon.nx_: the number of mass nodes along the width
% * _mmCon.ny_: the number of mass nodes along the height
% * _mmCon.n_: the total number of mass nodes
% * _mmCon.d_: the relative smoothing lenght of the kernel function
% associated to the mass distribution
% * _mmCon.rhoMin_: the minimum density (so that the stiffness matrix is
% not singular)
% * _mmCon.rhoMax_: the maximum density allowed in the asymptotic density
% * _mmCon.distrType_: the distribution pattern (1: regular, 2: random)
%
%
% The mass nodes can therefore be distributed in a rectangle of dimensions
% _mmCon.Lx * mmCon.Ly_ or randomly throughout the domain. The nodal masses
% _mmCon.mi_ are set so that the density is approximatively equal to 1
% inside the rectangle, and so that the total mass of the structure is
% approximatively equal to _mmCon.Lx * mmCon.Ly_.
%
% The function outputs are :
%
% * _mmCon_: a structure containing the mass distribution constants
% * _mnodes_: a cell structure containing the mass nodes coordinates,
% angles and dimensions
function [mmCon,mnodes] = massConstants(pCon,mCon)

    %Meshless mass constants
    mmCon.nx=6*pCon.Lx;                         % Number of mass nodes along the width
    mmCon.ny=3*pCon.Ly;                         % Number of mass nodes along the height
    mmCon.n=mmCon.nx*mmCon.ny;                  % Total number of mass nodes
    mmCon.d=1.5;                                % Relative smoothing length
    mmCon.m = mCon.m;                           % Number of integration cells
    mmCon.rhoMin = 1e-6;                        % Minimum density
    mmCon.rhoMax = 1.1;                         % Maximum density
    mmCon.distrType = 2;                        % Distribution type (1: in a rectangle, 2: random)
    
    
    % Nodes distribution parameters
    mmCon.Lx = pCon.Lx/2;                         % Rectangle length
    mmCon.Ly = pCon.Ly/6;                       % Rectangle height
    mmCon.drn = 0;%0.001*sqrt(mmCon.Lx^2 + mmCon.Ly^2);
    
    if mmCon.nx ~= 1
        mmCon.x0 = 0;                               % Rectangle low left corner x coordinate
        mmCon.dx=mmCon.Lx/(mmCon.nx-1);              % Horizontal distance between nodes
    else
        mmCon.x0 = pCon.Lx/2;
        mmCon.dx = pCon.Lx;
    end
    if mmCon.ny ~= 1
        mmCon.y0 = -mmCon.Ly/2;                      % Rectangle low left corner y coordinate
        mmCon.dy=mmCon.Ly/(mmCon.ny-1);              % Vertical distance between
    else
        mmCon.y0 = 0;
        mmCon.dy = pCon.Ly;
    end
    mmCon.dm=[mmCon.d*mmCon.dx ; mmCon.d*mmCon.dy];  % Smoothing length in x and y direction, respectively.
    mmCon.mi = mmCon.dx*mmCon.dy;               % Mass per node     

    % Create mass nodes
    mnodes = struct;
    for i=1:mmCon.n
        mnodes(i).x=[];
        mnodes(i).theta=0;
        mnodes(i).l = 2*[mmCon.dm(1);mmCon.dm(2)];
        mnodes(i).m = mmCon.mi;
    end
    
    % Mass nodes coordinates
    if mmCon.distrType == 1             % Regular distribution
        for i=1:mmCon.nx
            for j=1:mmCon.ny
                theta = rand()*2*pi;
                mnodes((i-1)*mmCon.ny+j).x=[mmCon.dx*(i-1)+mmCon.x0;...
                    mmCon.dy*(j-1)+mmCon.y0]+...
                    rand()*mmCon.drn()*[cos(theta);sin(theta)];
            end
        end
    elseif mmCon.distrType == 2         % Random distribution
        for i=1:mmCon.nx
            for j=1:mmCon.ny
                mnodes((i-1)*mmCon.ny+j).x= [rand()*pCon.Lx;...
                    rand()*pCon.Ly-pCon.Ly/2];
                mnodes((i-1)*mmCon.ny+j).theta = rand()*pi;
            end
        end
    end
    
end