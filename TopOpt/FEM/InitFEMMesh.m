%% Finite Element Method (FEM) mesh initialization
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

%% Discretization constants definition 

GlobalConst
 
tic %Mesh timer


mCon.mx=20*pCon.Lx;                         % Number of elements along the width
mCon.my=20*pCon.Ly;                         % Number of elements along the height
mCon.pn=1;                                  % Element degree (1: Q4, 2: Q8, 3: Q12)
mCon.nG=2;                                  % Number of quadrature points in one dimension in each cells

mCon.nx=mCon.pn*mCon.mx+1;                  % Number of discretization nodes along the width
mCon.ny=mCon.pn*mCon.my+1;                  % Number of discretization nodes along the height
mCon.nb=mCon.ny-1;                          % Number of boundary nodes per boundary
mCon.dx=pCon.Lx/(mCon.nx-1);                % Horizontal distance between nodes
mCon.dy=pCon.Ly/(mCon.ny-1);                % Vertical distance between nodes
mCon.cdx=pCon.Lx/mCon.mx;                   % Horizontal distance between cells
mCon.cdy=pCon.Ly/mCon.my;                   % Vertical distance between cells
mCon.n=4*mCon.pn + (3*mCon.pn-1)*...
    (mCon.mx+mCon.my-2) + (1+2*(mCon.pn-1))*...
    (mCon.mx-1)*(mCon.my-1);                % Total number of nodes
mCon.xc0 = 0;                               % Mesh x0 coord.
mCon.yc0 = -pCon.Ly/2;                      % Mesh y0 coord.

mCon.Lx = pCon.Lx;                          % Mesh length
mCon.Ly = pCon.Ly;                          % Mesh width
mCon.x0 = mCon.xc0;                         % Background mesh x0 coord
mCon.y0 = mCon.yc0;                         % Background mesh y0 coord
mCon.m=mCon.mx*mCon.my;                     % Number of total internal cells
mCon.bcdy=pCon.Ly/mCon.nb;                  % Vertical distance between nodes 
mCon.mbl=mCon.nb*pCon.nlLoad;               % Number of boundary cells along the right edge
mCon.mb=mCon.mbl;                           % Total number of line boundary cells
mCon.mp = pCon.npLoad;                      % Total number of point boundary cells



%% Nodes and cells creation

%Create nodes
nodes = struct;
for i = 1 : mCon.n
    nodes(i).x = [];
end

% Create elements
cells = struct;
for i=1:mCon.m
    cells(i).x=[];                  % Coordinates of elements center
    cells(i).dx=[];                 % Width and height of elements
    cells(i).ni=mCon.nG*mCon.nG;    % Number of integration points
    cells(i).J=[];                  % Jacobian
    cells(i).nen = [];              % Nodes of the element
    for j=1:cells(i).ni  
        cells(i).int(j).x=[];       % Coordinates of integration points 
        cells(i).int(j).nemn=[];    % Neighboring mass nodes
        cells(i).int(j).nefmn=[];   % Neighboring fixed mass nodes
        cells(i).int(j).w=[];       % Weight values for each intergration point
        cells(i).int(j).cv=[];      % Body force vector
    end
end

%Create boundary cells
bcells = struct;
for i=1:mCon.mb+mCon.mp
    bcells(i).BC=0;                 % Boundary type - 0 is undefined
    bcells(i).x=[];                 % Cells center
    bcells(i).dx=[];                % Width of cells 
    if i <= mCon.mb
        bcells(i).ni=mCon.nG;           % Number of integration points
    else
        bcells(i).ni=1;                 % Number of integration points
    end
    bcells(i).J=[];                 % Jacobian
    bcells(i).Ni=[];                % Two neighboring element numbers - only necessary for essential boundary condition
    bcells(i).n=[];                 % Two neighboring element numbers - only necessary for essential boundary condition
    for j=1:bcells(i).ni  
        bcells(i).int(j).x=[];         % Integration points coordinates
        bcells(i).int(j).nec=[];       % Containing cells
        bcells(i).int(j).w=[];         % Weight values for each intergration point
        bcells(i).int(j).bv=[];        % Boundary specific scalar or vector
    end
end

%% Distribution
k = 0;
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
nl = 1;
np = 1;
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
            bcells(i).int(k).x=pCon.lLoad(j).param(eta+t(k)*Gmx);
            bcells(i).int(k).w=w(k);
            tb = pCon.lLoad(j).F(eta+Gcx+Gmx*t(l));
            bcells(i).int(k).bv=tb;
        end
        if ~isempty(pCon.loads) && j > pCon.loads(1,nl)
            nl = nl + 1;
        end
        bcells(i).n=nl;
    elseif bcells(i).BC==4 
        j = i - pCon.nlLoad*mCon.nb;
        bcells(i).x = pCon.pLoad(j).x;
        bcells(i).J = 1;
        bcells(i).int(1).x = bcells(i).x;
        bcells(i).int(1).w = 1;
        bcells(i).int(1).bv = pCon.pLoad(j).F;
        if ~isempty(pCon.loads) && j > pCon.loads(2,np)
            np = np + 1;
        end
        bcells(i).n=np;
    end
end

%% Setting domains
[XX,YY] = meshgrid(linspace(0,pCon.Lx,mCon.mx+1),linspace(-pCon.Ly/2,pCon.Ly/2,mCon.my+1));
remove = [];
for d = 1: length(pCon.domains)
    f = pCon.domains(d).f(XX,YY);
    for i = 1 : mCon.mx
        for j = 1 : mCon.my
            if max(max(f(j:j+1,i:i+1))) < 0
                remove = [remove,(i-1)*mCon.my + mod(j-1,mCon.my)+1];
            end
        end
    end
end
    
%% Clearing holes    
    
for h = 1 : length(pCon.holes)
    if pCon.holes(h).type == 1                            % Rectangle
        xdHole = pCon.holes(h).x0 - pCon.holes(h).l/2;
        xuHole = pCon.holes(h).x0 + pCon.holes(h).l/2;
        for i = 1:mCon.m                            % Removing cells
            if pCon.holes(h).coverBoundary
                tmp = min(cells(i).x+cells(i).dx/2 <= xuHole) && min(cells(i).x-cells(i).dx/2 >= xdHole);
            else
                tmp = min(cells(i).x-cells(i).dx/2 < xuHole) && min(cells(i).x+cells(i).dx/2 > xdHole);
            end
            if tmp
                remove = [remove,i];
            end
        end
    elseif pCon.holes(h).type == 2
        for i = 1:mCon.m                            % Removing cells
            if pCon.holes(h).coverBoundary
                tmp = max([cells(i).x(1)-cells(i).dx(1)/2 ; cells(i).x(2)-cells(i).dx(2)/2] - pCon.holes(h).x0) <= pCon.holes(h).r && ...
                    max([cells(i).x(1)+cells(i).dx(1)/2 ; cells(i).x(2)-cells(i).dx(2)/2] - pCon.holes(h).x0) <= pCon.holes(h).r && ...
                    max([cells(i).x(1)-cells(i).dx(1)/2 ; cells(i).x(2)+cells(i).dx(2)/2] - pCon.holes(h).x0) <= pCon.holes(h).r && ...
                    max([cells(i).x(1)+cells(i).dx(1)/2 ; cells(i).x(2)+cells(i).dx(2)/2] - pCon.holes(h).x0) <= pCon.holes(h).r ;
            else
                tmp = norm([cells(i).x(1)-cells(i).dx(1)/2 ; cells(i).x(2)-cells(i).dx(2)/2] - pCon.holes(h).x0) <= pCon.holes(h).r || ...
                    norm([cells(i).x(1)+cells(i).dx(1)/2 ; cells(i).x(2)-cells(i).dx(2)/2] - pCon.holes(h).x0) <= pCon.holes(h).r || ...
                    norm([cells(i).x(1)-cells(i).dx(1)/2 ; cells(i).x(2)+cells(i).dx(2)/2] - pCon.holes(h).x0) <= pCon.holes(h).r || ...
                    norm([cells(i).x(1)+cells(i).dx(1)/2 ; cells(i).x(2)+cells(i).dx(2)/2] - pCon.holes(h).x0) <= pCon.holes(h).r ;
            end
            if tmp
                remove = [remove,i];
            end
        end
    end
end
remove = unique(remove);


% Removing nodes without cells

removeNodes = [];
neighboringCells = cell(mCon.n,1);
for n = 1 : mCon.n
   colInd = floor((n-1)/(mCon.my*(2*mCon.pn-1)+mCon.pn));
   nwoCol = mod(n-1,mCon.my*(2*mCon.pn-1)+mCon.pn);
   verticalNode = nwoCol <= mCon.pn*mCon.my;
   cornerNode = verticalNode && ~mod(nwoCol,mCon.pn);
   
   
   if cornerNode
       % Find at most four neighboring cells
       lineInd = nwoCol/mCon.pn;
       
       if colInd == 0
          
           if lineInd == 0
               cellInd = 1;
           elseif lineInd == mCon.my
               cellInd = mCon.my;
           else
               cellInd = lineInd + [0 ; 1];
           end
           
       elseif colInd == mCon.mx
           
           if lineInd == 0
               cellInd = mCon.m - mCon.my + 1;
           elseif lineInd == mCon.my
               cellInd = mCon.m;
           else
               cellInd = mCon.m - mCon.my + lineInd + [ 0 ; 1];
           end
           
       else
          
           if lineInd == 0
               cellInd = [1+(colInd-1)*mCon.my ; 1+colInd*mCon.my];
               
           elseif lineInd == mCon.my
               cellInd = [colInd*mCon.my ; (colInd+1)*mCon.my];
               
           else
               cellInd = [lineInd+(colInd-1)*mCon.my ; 1+lineInd+(colInd-1)*mCon.my ; ...
                   lineInd+colInd*mCon.my ; 1+lineInd+colInd*mCon.my];
           end
           
       end
   else
       % Find at most two neighboring cells
       if verticalNode
           lineInd = floor(nwoCol/mCon.pn);
            if colInd == 0 
                cellInd = 1+lineInd;
            elseif colInd == mCon.mx
                cellInd = 1+(mCon.mx-1)*mCon.my+lineInd;
            else
                cellInd = [1+(colInd-1)*mCon.my+lineInd ; 1+colInd*mCon.my+lineInd];
            end
       else
           lineInd = mod(nwoCol-(mCon.pn*mCon.my)-1,mCon.my+1);
           if lineInd == 0 
                cellInd = 1+colInd*mCon.my;
            elseif lineInd == mCon.my
                cellInd = colInd*mCon.my+lineInd;
            else
                cellInd = [colInd*mCon.my+lineInd ; 1+colInd*mCon.my+lineInd];
            end
       end
   end
   if ismember(cellInd,remove)
       removeNodes = [removeNodes,n];
   end
   neighboringCells{n} = cellInd;
end 
cells(remove) = [];
correspCellInd = zeros(1,mCon.m);
correspCellInd(setdiff(1:mCon.m,remove)) = 1:length(cells);
mCon.m = length(cells);



nodes(removeNodes) = [];
correspNodeInd = setdiff(1:mCon.n,removeNodes);
mCon.n = mCon.n - length(removeNodes);



% Volume computation
pCon.vol = mCon.m*prod(cells(1).dx);
   
%% Neighboring relations    
% Find element nodes
for i = 1 : mCon.n                        % i: new node index
    n = correspNodeInd(i);                                  % n: old node index
    for j = 1 : length(neighboringCells{n})
        nc = correspCellInd(neighboringCells{n}(j));        % nc: new cell index
        if nc ~= 0
            cells(nc).nen = [cells(nc).nen,i];
        end
    end
end

% Essential boundary conditions
tol = 1e-9;
x1=[nodes.x];
nodelabel=1:mCon.n;
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
       ubari = pCon.lbc(i).u(pCon.lbc(i).param(l));
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
x2 = [cells.x];
dx2 = [cells.dx]/2;
for j=1:mCon.mb+mCon.mp
    for k=1:bcells(j).ni
        bcells(j).int(k).nec=celllabel(and(and(bcells(j).int(k).x(1)>=x2(1,:)-dx2(1,:)-tol,bcells(j).int(k).x(1)<=x2(1,:)+dx2(1,:)+tol),...
            and(bcells(j).int(k).x(2)>=x2(2,:)-dx2(2,:)-tol,bcells(j).int(k).x(2)<= x2(2,:)+dx2(2,:)+tol)));
    end
end
    
ndofsE = 2*length(cells(1).nen);
lInd = mCon.m*ndofsE^2;

mCon.iK = zeros(lInd,1);
mCon.jK = zeros(lInd,1);
mCon.w = zeros(mCon.m*mCon.nG^2,1);

for ic = 1 : mCon.m
    en=zeros(ndofsE,1);
    en(1:2:end-1)=2*[cells(ic).nen]-1;
    en(2:2:end)=2*[cells(ic).nen];
    mCon.iK((ic-1)*ndofsE^2+1:ic*ndofsE^2) = repmat(en,ndofsE,1);
    mCon.jK((ic-1)*ndofsE^2+1:ic*ndofsE^2) = kron(en,ones(ndofsE,1));
    for ip=1:cells(ic).ni    
        mCon.w((ic-1)*mCon.nG^2+ip) = cells(ic).J*cells(ic).int(ip).w;
    end
end


disp('Mesh initialized')
disp(['     ',num2str(mCon.m),' elements    (',num2str(mCon.mx),' along x and ',num2str(mCon.my),' along y)'])
disp(['     ',num2str(mCon.n),' nodes       (',num2str(mCon.nx),' along x and ',num2str(mCon.ny),' along y)'])

%% Mass nodes creation

[mmCon,mnodes] = massConstants(pCon,mCon);

    
time1=toc; %Mesh timer
%disp([num2str(time1),' seconds to create problem'])

end
