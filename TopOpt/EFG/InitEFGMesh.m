%% Element-Free Galerkin (EFG) mesh initialization
% Created by J.T.B. Overvelde - 18 April 2012
%
% Master's thesis - The Moving Noda Approach in Topology Optimization
%
% <http://www.overvelde.com>
%
% Modified by G. Raze - June 2016
% 
% Initialize problem, nodal distribution and background mesh for an EFG
% analysis. The global constants are modified.

function time1=InitEFGMesh()

%% Constants definition 

GlobalConst
 
tic %Mesh timer

pCon = problemConstants();

%Meshless discretization constants
mCon.addCells = 0;                          % Background mesh domain 
                                            % multiplication factor
                                            % /!\ singularities if
                                            % mCon.mass = false
mCon.nx=10*pCon.Lx+1+2*mCon.addCells;                                 % Number of discretization nodes along the width
mCon.ny=10*pCon.Ly+1+2*mCon.addCells;                                 % Number of discretization nodes along the height
mCon.nb=mCon.ny-1;                          % Number of boundary nodes per boundary
mCon.dx=pCon.Lx/(mCon.nx-1);                %Horizontal distance between nodes
mCon.dy=pCon.Ly/(mCon.ny-1);                %Vertical distance between nodes
mCon.n=mCon.nx*mCon.ny;                     % Total number of nodes
mCon.xc0 = -mCon.addCells*mCon.dx;                    % Background mesh x0 coord.
mCon.yc0 = -pCon.Ly/2-mCon.addCells*mCon.dy;          % Background mesh y0 coord.

mCon.Lx = pCon.Lx+2*mCon.addCells*mCon.dx;           % Background mesh length
mCon.Ly = pCon.Ly+2*mCon.addCells*mCon.dy;           % Background mesh width
mCon.x0 = mCon.xc0;                                  % Background mesh x0 coord
mCon.y0 = mCon.yc0;                                  % Background mesh y0 coord
mCon.d=2.5;                                 %Relative smoothing length
mCon.pn=2;                                  %Number of monomials in the Monomial basis vector
mCon.nG=2;                                  %Number of quadrature points in one dimension in each cells
mCon.BCuType=1;                             %Displacement Boundary type, 1: integration cells 2: Collocated method
mCon.DispNodesType=1;                       %Nodal distribution type, 1: pattern 2: random 3: semi-random
mCon.mx=mCon.nx-1;                          %Number of cells along the width
mCon.my=mCon.ny-1;                          %Number of cells along the height
mCon.m=mCon.mx*mCon.my;                     % Number of total internal cells
mCon.bcdy=pCon.Ly/mCon.nb;                  % Vertical distance between nodes

if mCon.BCuType==1       
    mCon.mbu=mCon.nb*pCon.nlbc;             % Number of boundary cells along the left edge
    mCon.mbuc=0;                            % Number of collocated boundary particles along left edge
    mCon.me=mCon.mbu+pCon.nlbc+pCon.npbc;             % Number of essential boundary condition nodes
else
    mCon.mbu=0;                             % Number of boundary cells along the left edge
    mCon.mbuc=mCon.nb*pCon.nlbc;            % Number of collocated boundary particles along left edge
    mCon.me=0;                              % Number of essential boundary condition nodes
end
mCon.mbl=mCon.nb*pCon.nlLoad;              % Number of boundary cells along the right edge
mCon.mb=mCon.mbl+mCon.mbu;                  % Total number of line boundary cells
mCon.dm=[mCon.d*mCon.dx ; mCon.d*mCon.dy];  % Smoothing length in x and y direction, respectively.
mCon.mp = pCon.npbc+pCon.npLoad;            % Total number of point boundary cells
mCon.drn=1;                                 % Diameter of semi-random distribution

[mmCon,mnodes] = massConstants(pCon,mCon);

% %Position in domain for which the final solution is calculated
% sCon.snx=25;                                %Number of nodes in x direction for which the solution is calculated
% sCon.sny=25;                                %Number of nodes in y direction for which the solution is calculated
% sCon.sn=sCon.snx*sCon.sny;                  %Number of nodes for which the solution is calculated
% sCon.dsx=pCon.Lx/(sCon.snx-1);              %Horizontal distance between solution nodes
% sCon.dsy=pCon.Ly/(sCon.sny-1);              %Vertical distance between solution nodes



%% Nodes and cells creation

%Create nodes
nodes = struct;
for i=1:mCon.n
    nodes(i).x=[];       %nodes coordinates
    nodes(i).nen=[];     %Neighboring nodes
end
% %Create Solution nodes
% solnodes = struct;
% for i=1:sCon.sn
%     solnodes(i).x=[];    %solnodes coordinates
%     solnodes(i).nen=[];  %Neighboring nodes
%     solnodes(i).nx=[];   %Matrix x location
%     solnodes(i).ny=[];   %Matrix y location
% end


%Create internal background cells
cells = struct;
for i=1:mCon.m
    cells(i).x=[];                  %coordinates of cells center
    cells(i).dx=[];                 %width and height cells center
    cells(i).ni=mCon.nG*mCon.nG;    %number of integration points
    cells(i).J=[];                  %Jacobian
    for j=1:cells(i).ni  
        cells(i).int(j).x=[];       %Coordinates of integration points 
        cells(i).int(j).nen=[];     %Neighboring nodes
        cells(i).int(j).nemn=[];    %Neighboring mass nodes
        cells(i).int(j).w=[];       %Weight values for each intergration point
        cells(i).int(j).cv=[];      %Body force vector
    end
end
%Create boundary cells
bcells = struct;
for i=1:mCon.mb+mCon.mp
    bcells(i).BC=0;                 %Boundary type - 0 is undefined
    bcells(i).x=[];                 %cells center
    bcells(i).dx=[];                %width and height cells center
    if i <= mCon.mb
        bcells(i).ni=mCon.nG;           %number of integration points
    else
        bcells(i).ni=1;                 %number of integration points
    end
    bcells(i).J=[];                 %Jacobian
    bcells(i).Ni=[];                %Two neighboring element numbers - only necessary for essential boundary condition
    for j=1:bcells(i).ni  
        bcells(i).int(j).x=[];           %Integration points
        bcells(i).int(j).nen=[];         %Neighboring nodes
        bcells(i).int(j).nemn=[];        %Neighboring mass nodes
        bcells(i).int(j).w=[];           %Weight values for each intergration point
        bcells(i).int(j).bv=[];          %Boundary specific scalar or vector
        bcells(i).int(j).N=zeros(2,1);   %Interpolation value - only necessary for essential boundary condition
    end
end

%% Distribution

ag=1; %Constant telling if the material distribution is homogeneous
while ag==1
    ag=0;
    if mCon.DispNodesType==1 %Organized distribution of particles
        for i=1:mCon.nx
            for j=1:mCon.ny
                nodes((i-1)*mCon.ny+j).x=[mCon.dx*(i-1)+mCon.x0; mCon.dy*(j-1)+mCon.y0];
            end
        end
    elseif mCon.DispNodesType==2 %Random distribution in domain
        for i=1:mCon.n
            nodes(i).x=[mCon.Lx*rand(1)+mCon.x0; mCon.Ly*rand(1)+mCon.y0];
        end
    elseif mCon.DispNodesType==3 %First Organize the nodes, then distribute them randomly with small variation around node
        for i=1:mCon.nx
            for j=1:mCon.ny
                nodes((i-1)*mCon.ny+j).x=[mCon.dx*(i-1)+mCon.x0; mCon.dy*(j-1)+mCon.y0];
            end
        end
        for i=1:mCon.n
            a=1;
            while a==1
                a=0;
                dx=mCon.drn*[mCon.dx*(rand(1)-0.5); mCon.dy*(rand(1)-0.5)];
                if nodes(i).x(2)+dx(2)>mCon.Ly/2
                    a=1;
                elseif nodes(i).x(2)+dx(2)<-mCon.Ly/2
                    a=1;
                else
                    nodes(i).x=nodes(i).x+dx;
                end
            end
        end
    end
    


%     %solnodes coordinates
%     for i=1:sCon.snx
%         for j=1:sCon.sny
%             solnodes((i-1)*sCon.sny+j).x=[sCon.dsx*(i-1); sCon.dsy*(j-1)-pCon.Ly/2];
%             solnodes((i-1)*sCon.sny+j).nx=i;
%             solnodes((i-1)*sCon.sny+j).ny=j;
%         end
%     end

    %Boundary type
    for i=1:mCon.mb+mCon.mp
        if i<=mCon.mbu 
            bcells(i).BC=1; % Line essential boundary condition
        elseif i<=mCon.mb
            bcells(i).BC=2; % Line traction boundary condition
        elseif i<=mCon.mb+pCon.npbc
            bcells(i).BC=3; % Point essential boundary condition
        else
            bcells(i).BC=4; % Point traction boundary condition
        end
    end

    [t,w]=ConGauss(mCon.nG); %Gaus quadrature in general coordinates

    %Internal cells coordinates, integration points and weights
    for i=1:mCon.mx
        for j=1:mCon.my
            cells((i-1)*mCon.my+j).x=[mCon.dx/2+mCon.dx*(i-1)+mCon.xc0; mCon.dy/2+mCon.dy*(j-1)+mCon.yc0];
            cells((i-1)*mCon.my+j).dx=[mCon.dx;mCon.dy];   
            x1=cells((i-1)*mCon.my+j).x(1)-cells((i-1)*mCon.my+j).dx(1)/2;
            x2=cells((i-1)*mCon.my+j).x(1)+cells((i-1)*mCon.my+j).dx(1)/2;
            y1=cells((i-1)*mCon.my+j).x(2)-cells((i-1)*mCon.my+j).dx(2)/2;
            y2=cells((i-1)*mCon.my+j).x(2)+cells((i-1)*mCon.my+j).dx(2)/2;
            Gcx=1/2*(x2+x1); Gmx=1/2*(x2-x1);
            Gcy=1/2*(y2+y1); Gmy=1/2*(y2-y1);
            cells((i-1)*mCon.my+j).J=Gmx*Gmy;
            for k=1:mCon.nG
                for l=1:mCon.nG
                    cells((i-1)*mCon.my+j).int((k-1)*mCon.nG+l).x=[Gcx+Gmx*t(k);Gcy+Gmy*t(l)];
                    cells((i-1)*mCon.my+j).int((k-1)*mCon.nG+l).w=w(k)*w(l);
                    cells((i-1)*mCon.my+j).int((k-1)*mCon.nG+l).cv=pCon.b;
                end
            end
        end
    end
    
    %Boundary cells coordinates, integration points and weights
    for i=1:mCon.mb+mCon.mp
        if bcells(i).BC==1
            j = ceil(i/mCon.nb);
            k = mod(i-1,mCon.nb);
            dx = pCon.lbc(j).length/mCon.nb;
            eta = k*dx+dx/2;
            bcells(i).x=pCon.lbc(j).param(eta);
            bcells(i).dx=dx;
            bcells(i).Ni=[i+(j-1) i+(j-1)+1];   
            x1=-bcells(i).dx/2;
            x2=+bcells(i).dx/2;
            Gcx=1/2*(x2+x1); Gmx=1/2*(x2-x1);
            bcells(i).J=Gmx;
            for l=1:mCon.nG
                bcells(i).int(l).x=bcells(i).x+Gcx+Gmx*t(l)*pCon.lbc(j).dir;
                bcells(i).int(l).w=w(l);
                bcells(i).int(l).N=[1-(t(l)+1)/2, (t(l)+1)/2];
                ub = pCon.lbc(j).u(eta+Gcx+Gmx*t(l));
                bcells(i).int(l).bv=ub;
            end
        elseif bcells(i).BC==2
            j = ceil(i/mCon.nb)-pCon.nlbc;
            k = mod(i-1,mCon.nb);
            dx = pCon.lLoad(j).length/mCon.nb;
            eta = k*dx+dx/2;
            bcells(i).x=pCon.lLoad(j).param(eta);
            bcells(i).dx=dx;
            x1=-bcells(i).dx/2;
            x2=+bcells(i).dx/2;
            Gcx=1/2*(x2+x1); Gmx=1/2*(x2-x1);
            bcells(i).J=Gmx;
            for k=1:mCon.nG
                bcells(i).int(k).x=bcells(i).x+[0;Gcx+Gmx*t(k)];
                bcells(i).int(k).w=w(k);
                tb = pCon.lLoad(j).F(eta+Gcx+Gmx*t(l));
                bcells(i).int(k).bv=tb;
            end
        elseif bcells(i).BC==3    
            j = i - (pCon.nlbc+pCon.nlLoad)*mCon.nb;
            bcells(i).x = pCon.pbc(j).x;
            bcells(i).Ni = [i-pCon.nlLoad*(mCon.nb-1) i-pCon.nlLoad*(mCon.nb-1)+1];
            bcells(i).J = 1;
            bcells(i).int(1).x = bcells(i).x;
            bcells(i).int(1).w = 2;
            bcells(i).int(1).N=[1/2, 1/2];
            bcells(i).int(1).bv = pCon.pbc(j).u;
        elseif bcells(i).BC==4 
            j = i - (pCon.nlbc+pCon.nlLoad)*mCon.nb - pCon.npbc;
            bcells(i).x = pCon.pLoad(j).x;
            bcells(i).J = 1;
            bcells(i).int(1).x = bcells(i).x;
            bcells(i).int(1).w = 2;
            bcells(i).int(1).bv = pCon.pLoad(j).F;
        end
    end
    bcCon.Bcux=[]; bcCon.Bcuy=[]; bcCon.BCunx=0; bcCon.BCuny=0;
    if mCon.BCuType~=1 %Collocated BC for diplacement
        rep=0;
        for i=1:mCon.mblc
            rep=rep+1;
            u=AnalSol(nodes(i).x);
            bcCon.BCux(1,rep)=i;
            bcCon.BCux(2,rep)=u(1);
            bcCon.BCuy(1,rep)=i;
            bcCon.BCuy(2,rep)=u(2);
        end
        bcCon.BCunx=size(bcCon.BCux,2);
        bcCon.BCuny=size(bcCon.BCuy,2);
    end
%% Neighboring relations    

    
    x1=[nodes.x];
    nodelabel=1:mCon.n;
    %Find integraion points in influence domain of nodes
    for j=1:mCon.m
        for k=1:cells(j).ni
            cells(j).int(k).nen=nodelabel(and(abs(x1(1,:)-cells(j).int(k).x(1))<mCon.dm(1),abs(x1(2,:)-cells(j).int(k).x(2))<mCon.dm(2)));
%             if length(cells(j).int(k).nen)<5
%                 ag=1;
%             end    
        end
    end
    %Find boundary integration points in influence domain of nodes
    for j=1:mCon.mb+mCon.mp
        for k=1:bcells(j).ni
            bcells(j).int(k).nen=nodelabel(and(abs(x1(1,:)-bcells(j).int(k).x(1))<mCon.dm(1),abs(x1(2,:)-bcells(j).int(k).x(2))<mCon.dm(2)));
%             if length(bcells(j).int(k).nen)<5
%                 ag=1;
%             end    
        end
    end
    
    %Find neighbours of nodes
    for i=1:mCon.n
        nodes(i).nen=nodelabel(and(abs(x1(1,:)-nodes(i).x(1))<mCon.dm(1),abs(x1(2,:)-nodes(i).x(2))<mCon.dm(2)));
    end
    
%     %Find neighbour nodes of solnodes
%     for i=1:sCon.sn
%         solnodes(i).nen=nodelabel(and(abs(x1(1,:)-solnodes(i).x(1))<mCon.dm(1),abs(x1(2,:)-solnodes(i).x(2))<mCon.dm(2)));
%     end
%     
end

time1=toc; %Mesh timer
%disp([num2str(time1),' seconds to create problem'])
end
