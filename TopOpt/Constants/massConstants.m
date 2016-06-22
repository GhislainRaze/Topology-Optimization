%% Mass Constants
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
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
% * _mmCon.EMin_: the minimum Young modulus (so that the stiffness matrix 
% is not singular)
% * _mmCon.rhoMax_: the maximum density allowed in the asymptotic density
% * _mmCon.distrType_: the distribution pattern (1: regular, 2: random)
%
%
% The mass nodes can therefore be distributed in a rectangle of dimensions
% _mmCon.Lx * mmCon.Ly_, randomly throughout the domain, or randomly
% shifted from their position in the rectangle by a maximum radius
% _mmCon.drn_. The nodal masses _mmCon.mi_ are set so that the density is
% approximately equal to 1 inside the rectangle, and so that the total mass
% of the structure is approximately equal to _mmCon.Lx * mmCon.Ly_. The
% volume of the structure is given by _mmCon.vol_.
%
% The function outputs are :
%
% * _mmCon_: a structure containing the mass distribution constants
% * _mnodes_: a cell structure containing the mass nodes coordinates,
% angles and dimensions

function [mmCon,mnodes] = massConstants(pCon,mCon)

    % Meshless mass constants
    mmCon.nx=3*pCon.Lx;                             % Number of mass nodes along the width
    mmCon.ny=2*pCon.Ly;                             % Number of mass nodes along the height
    mmCon.n=mmCon.nx*mmCon.ny;                      % Total number of mass nodes
    mmCon.d=1.5;                                    % Relative smoothing length
    mmCon.m = mCon.m;                               % Number of integration cells
    mmCon.EMin = 1e-9;                              % Minimum Young's modulus
    mmCon.rhoMax = 1.05;                             % Maximum density
    mmCon.distrType = 1;                            % Distribution type (1: in a rectangle, 2: random, 3: semi-random)
    

    
    % Nodes distribution parameters
    mmCon.Lx = 3*pCon.Lx/4;                             % Rectangle length
    mmCon.Ly = pCon.Ly/3;                           % Rectangle height
    mmCon.drn = 0.001*sqrt(mmCon.Lx^2 + mmCon.Ly^2);% Semi-random maximal radius
    mmCon.vol = mmCon.Lx*mmCon.Ly;                  % Volume of the structure
    
    if mmCon.nx ~= 1
        mmCon.x0 = (pCon.Lx-mmCon.Lx)/2;            % Rectangle low left corner x coordinate
        mmCon.dx=mmCon.Lx/(mmCon.nx-1);             % Horizontal distance between nodes
    else
        mmCon.x0 = pCon.Lx/2;
        mmCon.dx = pCon.Lx;
    end
    if mmCon.ny ~= 1
        mmCon.y0 = -mmCon.Ly/2;           % Rectangle low left corner y coordinate
        mmCon.dy=mmCon.Ly/(mmCon.ny-1);             % Vertical distance between
    else
        mmCon.y0 = 0;
        mmCon.dy = pCon.Ly;
    end
    
    mmCon.dm=[mmCon.d*mmCon.dx ; mmCon.d*mmCon.dy]; % Smoothing length in x and y direction, respectively.
    mmCon.mi = mmCon.dx*mmCon.dy;                   % Mass per node     
    mmCon.rm = 1/mmCon.d^2;                         % Ratio between nodal influence domain and mass
    mmCon.mMax = mmCon.n*mmCon.mi;                  % Total mass of the structure
    
    % Create mass nodes
    mnodes = struct;
    for i=1:mmCon.n
        mnodes(i).x=[];
        mnodes(i).theta=0;
        mnodes(i).l = 2*[mmCon.dm(1);mmCon.dm(2)];
        mnodes(i).m = mmCon.mi;
    end
    
    if ~isempty(pCon.holes)
        nodesInHoles = true;
    else
        nodesInHoles = false;
    end
    
    % Mass nodes coordinates
    if mmCon.distrType == 1             % Regular distribution
        for i=1:mmCon.nx
            for j=1:mmCon.ny
                mnodes((i-1)*mmCon.ny+j).x=[mmCon.dx*(i-1)+mmCon.x0;...
                    mmCon.dy*(j-1)+mmCon.y0];
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
    elseif mmCon.distrType == 3         % Semi-random distribution
        for i=1:mmCon.nx
            for j=1:mmCon.ny
                theta = rand()*2*pi;
                mnodes((i-1)*mmCon.ny+j).x=[mmCon.dx*(i-1)+mmCon.x0;...
                    mmCon.dy*(j-1)+mmCon.y0]+...
                    rand()*mmCon.drn()*[cos(theta);sin(theta)];
            end
        end
    elseif mmCon.distrType == 4
        dx = pCon.Lx/mmCon.nx;
        dy = pCon.Ly/mmCon.ny;
        posx = 0:mmCon.nx-1;
        posy = 0:mmCon.ny-1;
        [x,y] = meshgrid(dx/2+dx*posx,-pCon.Ly/2+dy/2+dy*posy);
        nxtmp = mmCon.nx;
        nytmp = mmCon.ny;
        incry = true;
        nbNodes = 0;
        while nodesInHoles
            nodesInHoles = false;
            nbNodes = 0;
        	for i = 1 : size(x,1)
            	for j = 1 : size(x,2)
                    if checkHolesNode([x(i,j);y(i,j)],pCon.holes);
                        nbNodes = nbNodes + 1;
                    end
                end
            end
            if size(x,1)*size(x,2)-nbNodes < mmCon.n
            	nodesInHoles = true;
                if incry
                    posx = 0:nxtmp;
                	nxtmp = nxtmp+1;
                    dx = pCon.Lx/nxtmp;
                    incry = false;
                else
                    posy = 0:nytmp;
                    nytmp = nytmp+1;
                    dy = pCon.Ly/nytmp;
                    incry = true;
                end
                [x,y] = meshgrid(dx/2+dx*posx,-pCon.Ly/2+dy/2+dy*posy);
            end    
        end
        
        k = 0;
        for j = 1 : size(x,2)
            for i = 1 : size(x,1)
                if ~checkHolesNode([x(i,j);y(i,j)],pCon.holes)
                    k = k +1;
                    if k > mmCon.n
                        break
                    end
                    mnodes(k).x = [x(i,j);y(i,j)];
                end
            end
            if k > mmCon.n
                break
            end
        end
    end
    
    
    % Check that there are no mass nodes in holes
    
    nodesInd = 1:mmCon.n;
    nodesMoved = [];
    
    while nodesInHoles && mmCon.distrType < 4
        for i = 1:length(nodesInd)
            nodesInHoles = false;
            if checkHolesNode(mnodes(nodesInd(i)).x,pCon.holes)
                mnodes(nodesInd(i)).x = [rand()*pCon.Lx;...
                    rand()*pCon.Ly-pCon.Ly/2];
                nodesInHoles = true;
                nodesMoved = [nodesMoved,i];
            end
        end
        nodesInd = nodesInd(nodesMoved);
        nodesMoved = [];
    end
    
    
end