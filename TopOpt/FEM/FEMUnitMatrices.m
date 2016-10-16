%% Finite Element Method (FEM): unit matrices
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Determines the stiffness matrices associated with a unit density at the 
% Gauss points for a given nodal distribution and background mesh. Also
% determines the nodal force vector and the imposed nodal displacements.
%
% The outputs are
%
% * _Ke_: the unit stiffness matrices (a cell of sparse matrices)
% * _f_: the nodal force vector
% * _ubar_: the imposed nodal displacements (set to NaN if there is no
% imposed displacement associated to a given node)
% * _time1_: the time to assemble the problem

function [Ke,f,ubar,time1]=FEMUnitMatrices()

GlobalConst
InitFEMMesh;
tic %Assembly timer

% Assembly K matrix and f vector internal cells
Ke = cell(mCon.nG^2,1);
if pCon.type == 1
    f=zeros(2*mCon.n,1);
elseif pCon.type == 2
    f=zeros(2*mCon.n,2);
end
Fe=zeros(mCon.nG^2,2*length(cells(1).nen));

% Compute the element unit stiffness matrix
for ip=1:cells(1).ni                               % Iterations over the cell Gauss points
    B=zeros(3,2*length(cells(1).nen));               
    coord = 2*(cells(1).int(ip).x-cells(1).x)./cells(1).dx;
    [phi,dphidx,dphidy]=FEMShape(coord,length(cells(1).nen));
    B(1,1:2:end-1)=2/cells(1).dx(1)*dphidx;
    B(2,2:2:end)=2/cells(1).dx(2)*dphidy;
    B(3,1:2:end-1)=2/cells(1).dx(2)*dphidy;
    B(3,2:2:end)=2/cells(1).dx(1)*dphidx;
    Fe(ip,1:2:end-1) = phi*cells(1).int(ip).w*cells(1).J;       
    Fe(ip,2:2:end) = Fe(ip,1:2:end-1); 
    Ke{ip} = B'*pCon.D*B*cells(1).int(ip).w*cells(1).J;
end

for ic=1:mCon.m                                         % Iterations over the internal cells
    en=zeros(1,2*length(cells(ic).nen));
    en(1:2:end-1)=2*[cells(ic).nen]-1; 
    en(2:2:end)=2*[cells(ic).nen];
    for ip=1:cells(ic).ni                               % Body forces
        f(en(1:2:end-1),:) =f(en(1:2:end-1),:)+kron(Fe(ip,1:2:end-1)'*cells(ic).int(ip).cv(1),ones(1,size(f,2)));
        f(en(2:2:end),:)=f(en(2:2:end),:)+kron(Fe(ip,2:2:end)'*cells(ic).int(ip).cv(2),ones(1,size(f,2)));
    end
end


%Assembly f vector
for ic=1:mCon.mb+mCon.mp
    if bcells(ic).BC==2 || bcells(ic).BC==4
        for ip=1:bcells(ic).ni
            nc = bcells(ic).int(ip).nec(1);
            T=zeros(2*length(cells(nc).nen),1);
            en=zeros(1,2*length(cells(nc).nen));
            coord = 2*(bcells(ic).int(ip).x-cells(nc).x)./cells(nc).dx;
            phi=FEMShape(coord,length(cells(nc).nen));
            for neni=1:length(cells(nc).nen)
                T(2*neni-1:2*neni)=phi(neni).*bcells(ic).int(ip).bv;
                en(2*neni-1:2*neni)=[2*cells(nc).nen(neni)-1; 2*cells(nc).nen(neni)];
            end
            f(en,bcells(ic).n)=f(en,bcells(ic).n)+T*bcells(ic).int(ip).w*bcells(ic).J;
        end  
    end
end

% Imposed displacement
ubar = nan(2*mCon.n,1);
for i = 1 : size(bnodes,2)
    
    indx = 2*bnodes(1,i)-1;
    indy = 2*bnodes(1,i);
    if ~isnan(bnodes(2,i))
        ubar(indx) = bnodes(2,i);
    end
    if ~isnan(bnodes(3,i))
        ubar(indy) = bnodes(3,i);
    end
end
time1=toc; %Assembly timer

end