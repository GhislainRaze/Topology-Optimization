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

%Assembly K matrix and f vector internal cells
Ke = cell(mCon.m*mCon.nG^2,1);
f=zeros(2*mCon.n,1);

for ic=1:mCon.m                                         % Iterations over the internal cells
    for ip=1:cells(ic).ni                               % Iterations over the cell Gauss points
        B=zeros(3,2*length(cells(ic).nen));             
        F=zeros(2*length(cells(ic).nen),1);
        en=zeros(1,2*length(cells(ic).nen));
        coord = 2*(cells(ic).int(ip).x-cells(ic).x)./cells(ic).dx;
        [phi,dphidx,dphidy]=FEMShape(coord,length(cells(ic).nen));
        B(1,1:2:end-1)=2/cells(ic).dx(1)*dphidx;
        B(2,2:2:end)=2/cells(ic).dx(2)*dphidy;
        B(3,1:2:end-1)=2/cells(ic).dx(2)*dphidy;
        B(3,2:2:end)=2/cells(ic).dx(1)*dphidx;
        F(1:2:end-1)=phi*cells(ic).int(ip).cv(1);       
        F(2:2:end)=phi*cells(ic).int(ip).cv(2);
        en(1:2:end-1)=2*[cells(ic).nen]-1;              % x index of neighboring cells
        en(2:2:end)=2*[cells(ic).nen];                  % y index
        Ke{(ic-1)*mCon.nG^2+ip} = B'*pCon.D*B*cells(ic).int(ip).w*cells(ic).J;
        f(en)=f(en)+F*cells(ic).int(ip).w*cells(ic).J;
    end
end


%Assembly G matrix and q vector boundary cells
for ic=1:mCon.mb+mCon.mp
    if bcells(ic).BC==2 || bcells(ic).BC==4
        for ip=1:bcells(ic).ni
            ncc = bcells(ic).int(ip).nec;
            for i = 1 : length(ncc)
                nc = ncc(i);
                T=zeros(2*length(cells(nc).nen),1);
                en=zeros(1,2*length(cells(nc).nen));
                coord = 2*(bcells(ic).int(ip).x-cells(nc).x)./cells(nc).dx;
                phi=FEMShape(coord,length(cells(nc).nen));
                for neni=1:length(cells(nc).nen)
                    T(2*neni-1:2*neni)=phi(neni).*bcells(ic).int(ip).bv;
                    en(2*neni-1:2*neni)=[2*cells(nc).nen(neni)-1; 2*cells(nc).nen(neni)];
                end
                f(en)=f(en)+T*bcells(ic).int(ip).w*bcells(ic).J;
            end
        end  
    end
end

% Imposed displacement
ubar = nan(2*mCon.n,1);
for i = 1 : length(bnodes)
    
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