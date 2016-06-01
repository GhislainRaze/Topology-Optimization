% Finite Element Method (FEM) mesh initialization
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Initialize problem and mesh for a FEM analysis. The global constants are
% modified.

function time1=InitFEMMesh()

%% Constants definition 

GlobalConst
 
tic %Mesh timer

pCon = problemConstants();

%Meshless discretization constants
mCon.addCells = 0;                          % Background mesh domain 
                                            % multiplication factor
mCon.mx=9*pCon.Lx;                          %Number of elements along the width
mCon.my=6*pCon.Ly;                          %Number of elements along the height
mCon.pn=1;                                  %Element degree (1: Q4, 2: Q8, 3: Q12)
mCon.nx=mCon.pn*mCon.mx+1+2*mCon.addCells;  % Number of discretization nodes along the width
mCon.ny=mCon.pn*mCon.my+1+2*mCon.addCells;  % Number of discretization nodes along the height
mCon.nb=mCon.ny-1;                          % Number of boundary nodes per boundary
mCon.dx=pCon.Lx/(mCon.nx-1);                %Horizontal distance between nodes
mCon.dy=pCon.Ly/(mCon.ny-1);                %Vertical distance between nodes
mCon.cdx=pCon.Lx/mCon.mx;                   %Horizontal distance between cells
mCon.cdy=pCon.Ly/mCon.my;                   %Vertical distance between cells
mCon.n=4*mCon.pn + (3*mCon.pn-1)*...
    (mCon.mx+mCon.my-2) + (1+2*(mCon.pn-1))*...
    (mCon.mx-1)*(mCon.my-1);                % Total number of nodes
mCon.xc0 = -mCon.addCells*mCon.dx;                    % Background mesh x0 coord.
mCon.yc0 = -pCon.Ly/2-mCon.addCells*mCon.dy;          % Background mesh y0 coord.

mCon.Lx = pCon.Lx+2*mCon.addCells*mCon.dx;           % Background mesh length
mCon.Ly = pCon.Ly+2*mCon.addCells*mCon.dy;           % Background mesh width
mCon.x0 = mCon.xc0;                                  % Background mesh x0 coord
mCon.y0 = mCon.yc0;                                  % Background mesh y0 coord
mCon.nG=2;                                  %Number of quadrature points in one dimension in each cells
mCon.DispNodesType=1;                       %Nodal distribution type, 1: pattern 2: random 3: semi-random
mCon.m=mCon.mx*mCon.my;                     % Number of total internal cells
mCon.bcdy=pCon.Ly/mCon.nb;                  % Vertical distance between nodes 
mCon.mbl=mCon.nb*pCon.nlLoad;               % Number of boundary cells along the right edge
mCon.mb=mCon.mbl;                           % Total number of line boundary cells
mCon.mp = pCon.npLoad;                      % Total number of point boundary cells
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
    cells(i).nen = [];              % Nodes of the element
    for j=1:cells(i).ni  
        cells(i).int(j).x=[];       %Coordinates of integration points 
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
        bcells(i).int(j).nec=[];       %Containing cells
        bcells(i).int(j).w=[];           %Weight values for each intergration point
        bcells(i).int(j).bv=[];          %Boundary specific scalar or vector
    end
end

%% Distribution

ag=1; %Constant telling if the material distribution is homogeneous
k = 0;
while ag==1
    ag=0;
    if mCon.DispNodesType==1 %Organized distribution of particles
        for i=1:mCon.nx
            for j=1:mCon.ny
                if ~mod(i-1,mCon.pn)
                    k = k + 1;
                    nodes(k).x=[mCon.dx*(i-1)+mCon.x0; mCon.dy*(j-1)+mCon.y0];
                elseif ~mod(j-1,mCon.pn)
                    k = k + 1;
                    nodes(k).x=[mCon.dx*(i-1)+mCon.x0; mCon.dy*(j-1)+mCon.y0];
                end
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
        if i<=mCon.mb
            bcells(i).BC=2; % Line traction boundary condition
        else
            bcells(i).BC=4; % Point traction boundary condition
        end
    end

    [t,w]=ConGauss(mCon.nG); %Gaus quadrature in general coordinates

    %Internal cells coordinates, integration points and weights
    for i=1:mCon.mx
        for j=1:mCon.my
            cells((i-1)*mCon.my+j).x=[mCon.cdx/2+mCon.cdx*(i-1)+mCon.xc0; mCon.cdy/2+mCon.cdy*(j-1)+mCon.yc0];
            cells((i-1)*mCon.my+j).dx=[mCon.cdx;mCon.cdy];   
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
        if bcells(i).BC==2
            j = ceil(i/mCon.nb);
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
        elseif bcells(i).BC==4 
            j = i - pCon.nlLoad*mCon.nb;
            bcells(i).x = pCon.pLoad(j).x;
            bcells(i).J = 1;
            bcells(i).int(1).x = bcells(i).x;
            bcells(i).int(1).w = 2;
            bcells(i).int(1).bv = pCon.pLoad(j).F;
        end
    end
   
%% Neighboring relations    
    
    tol = 1e-9;
    
    x1=[nodes.x];
    nodelabel=1:mCon.n;
    %Find element nodes
    for j=1:mCon.m
       cells(j).nen=nodelabel(and(and(x1(1,:)>=cells(j).x(1)-cells(j).dx(1)/2-tol,x1(1,:)<= cells(j).x(1)+cells(j).dx(1)/2+tol),...
            and(x1(2,:)>=cells(j).x(2)-cells(j).dx(2)/2-tol,x1(2,:)<= cells(j).x(2)+cells(j).dx(2)/2+tol)));
    end
    
    % Essential boundary conditions
    
    bnodes = [];
    
    for i=1:pCon.nlbc
       % Boundary along y
       newbnodes = nodelabel(and(and(pCon.lbc(i).x(1,1) == ...
           pCon.lbc(i).x(2,1),pCon.lbc(i).x(1,1) == x1(1,:)),...
           and(max(pCon.lbc(i).x(1,2),pCon.lbc(i).x(2,2)) >= x1(2,:),...
                min(pCon.lbc(i).x(1,2),pCon.lbc(i).x(2,2)) <= x1(2,:))));
       % Imposed displacement law
       for j = 1 : length(newbnodes)
           l = norm(nodes(newbnodes(j)).x' - pCon.lbc(i).x(1,:));
           ubari = pCon.lbc(i).u(l);
           bnodes = [bnodes, [newbnodes(j) ; ubari(1) ; ubari(2)]];
       end
       % Boundary along x
       newbnodes = nodelabel(and(and(pCon.lbc(i).x(1,2) == ...
           pCon.lbc(i).x(2,2),pCon.lbc(i).x(1,2) == x1(2,:)),...
           and(max(pCon.lbc(i).x(1,1),pCon.lbc(i).x(2,1)) >= x1(1,:),...
                min(pCon.lbc(i).x(1,1),pCon.lbc(i).x(2,1)) <= x1(1,:))));
            
       % Imposed displacement law
       for j = 1 : length(newbnodes)
           l = norm(nodes(newbnodes(j)).x' - pCon.lbc(i).x(1,:));
           ubari = pCon.lbc(i).param(l);
           bnodes = [bnodes, [newbnodes(j) ; ubari(1) ; ubari(2)]];
       end
    end
    
    
        
    
    for i = 1 : pCon.npbc
        newbnodes = nodelabel(and(pCon.pbc(i).x(1) == x1(1,:),...
            pCon.pbc(i).x(2) == x1(2,:)));
        for j = 1 : length(newbnodes)
           ubari = pCon.pbc(i).u;
           bnodes = [bnodes, [newbnodes(j) ; ubari(1) ; ubari(2)]];
       end
    end
    %Find boundary integration points in influence domain of nodes
    
    celllabel = 1:mCon.m;
    tcelllabel =  mCon.my:mCon.my:mCon.m;
    rcelllabel = mCon.m-mCon.my+1:mCon.m;
    x2 = [cells.x];
    dx2 = [cells.dx]/2;
    for j=1:mCon.mb+mCon.mp
        for k=1:bcells(j).ni
            bcells(j).int(k).nec=celllabel(and(and(bcells(j).int(k).x(1)>=x2(1,:)-dx2(1,:),bcells(j).int(k).x(1)<x2(1,:)+dx2(1,:)),...
                and(bcells(j).int(k).x(2)>=x2(2,:)-dx2(2,:),bcells(j).int(k).x(2)< x2(2,:)+dx2(2,:))));
            % Top boundary
            bcells(j).int(k).nec=[bcells(j).int(k).nec , tcelllabel(and(and(bcells(j).int(k).x(1)>=x2(1,mCon.my:mCon.my:end)-dx2(1,mCon.my:mCon.my:end),...
                bcells(j).int(k).x(1)<x2(1,mCon.my:mCon.my:end)+dx2(1,mCon.my:mCon.my:end)),...
                abs(bcells(j).int(k).x(2)-(x2(2,mCon.my:mCon.my:end)+dx2(2,mCon.my:mCon.my:end)))<tol))];
            % Right boundary
            bcells(j).int(k).nec=[bcells(j).int(k).nec , rcelllabel(and(and(bcells(j).int(k).x(2)>=x2(2,end-mCon.my+1:end)-dx2(2,end-mCon.my+1:end),...
                bcells(j).int(k).x(2)<x2(2,end-mCon.my+1:end)+dx2(2,end-mCon.my+1:end)),...
                abs(bcells(j).int(k).x(1)-(x2(1,end-mCon.my+1:end)+dx2(1,end-mCon.my+1:end)))<tol))];
            % Top right corner
            bcells(j).int(k).nec=[bcells(j).int(k).nec , mCon.m(and(abs(bcells(j).int(k).x(1)- x2(1,end)+dx2(1,end))<tol,...
                abs(bcells(j).int(k).x(2)- x2(2,end)+dx2(2,end))<tol))];
        end
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
