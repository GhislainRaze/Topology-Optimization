%% Element-Free Galerkin (EFG) method : Unit matrices
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
% determines the nodal force vector and the quantities associated to
% lagrangian multipliers.
%
% The outputs are
%
% * _Ke_: the unit stiffness matrices (a cell of sparse matrices)
% * _f_: the nodal force vector
% * _G_: the lagrangian multipliers matrix
% * _q_: the second member associated with lagrangian multipliers
% * _K_: the minimum stiffness matrix
% * _time1_: the time to assemble the problem

function [Ke,f,G,q,K,time1]=EFGUnitMatrices()

GlobalConst
InitEFGMesh;
tic %Assembly timer

%Assembly K matrix and f vector internal cells
Ke = cell(mCon.m*mCon.nG^2,1);
f=zeros(2*mCon.n,1);
K=zeros(2*mCon.n);

for ic=1:mCon.m                                             % Iterations over the internal cells
    for ip=1:cells(ic).ni                                   % Iterations over the cell Gauss points
        if ~isempty(cells(ic).int(ip).nen)
            B=zeros(3,2*length(cells(ic).int(ip).nen));     
            F=zeros(2*length(cells(ic).int(ip).nen),1);
            en=zeros(1,2*length(cells(ic).int(ip).nen));
            [phi,dphidx,dphidy]=MLSShape([nodes(cells(ic).int(ip).nen).x]',cells(ic).int(ip).x,mCon.dm,mCon.pn);
            B(1,1:2:end-1)=dphidx;
            B(2,2:2:end)=dphidy;
            B(3,1:2:end-1)=dphidy;
            B(3,2:2:end)=dphidx;
            F(1:2:end-1)=phi*cells(ic).int(ip).cv(1);       
            F(2:2:end)=phi*cells(ic).int(ip).cv(2);
            en(1:2:end-1)=2*[cells(ic).int(ip).nen]-1;      % x index of neighboring cells
            en(2:2:end)=2*[cells(ic).int(ip).nen];          % y index
            Ke{(ic-1)*mCon.nG^2+ip} = B'*pCon.D*B*cells(ic).int(ip).w*cells(ic).J;
            K(en,en) = K(en,en) + mmCon.EMin*B'*pCon.D*B*cells(ic).int(ip).w*cells(ic).J;
            f(en)=f(en)+F*cells(ic).int(ip).w*cells(ic).J;
        end
    end
end

K = sparse(K);

%Assembly G matrix and q vector boundary cells
G=zeros(2*mCon.n,2*mCon.me+bcCon.BCunx+bcCon.BCuny); 
q=zeros(2*mCon.me+bcCon.BCunx+bcCon.BCuny,1);
for ic=1:mCon.mb+mCon.mp
    if bcells(ic).BC==1 || bcells(ic).BC==3
        for in=1:2 %two-nodes element
            for ip=1:bcells(ic).ni
                U=zeros(2*length(bcells(ic).int(ip).nen),2);
                en=zeros(1,2*length(bcells(ic).int(ip).nen));
                [phi,dphidx,dphidy]=MLSShape([nodes(bcells(ic).int(ip).nen).x]',bcells(ic).int(ip).x,mCon.dm,mCon.pn);
                for neni=1:length(bcells(ic).int(ip).nen)
                    U(2*neni-1:2*neni,:)=phi(neni)*[bcells(ic).int(ip).N(in) 0; 0 bcells(ic).int(ip).N(in)];
                    en(2*neni-1:2*neni)=[2*bcells(ic).int(ip).nen(neni)-1; 2*bcells(ic).int(ip).nen(neni)];
                end
                q(2*bcells(ic).Ni(in)-1:2*bcells(ic).Ni(in))=q(2*bcells(ic).Ni(in)-1:2*bcells(ic).Ni(in))-[bcells(ic).int(ip).N(in) 0; 0 bcells(ic).int(ip).N(in)]*bcells(ic).int(ip).bv*bcells(ic).int(ip).w*bcells(ic).J;
                G(en,2*bcells(ic).Ni(in)-1:2*bcells(ic).Ni(in))=G(en,2*bcells(ic).Ni(in)-1:2*bcells(ic).Ni(in))-U*bcells(ic).int(ip).w*bcells(ic).J;
            end
        end
    elseif bcells(ic).BC==2 || bcells(ic).BC==4
        for ip=1:bcells(ic).ni
            T=zeros(2*length(bcells(ic).int(ip).nen),1);
            en=zeros(1,2*length(bcells(ic).int(ip).nen));
            [phi,dphidx,dphidy]=MLSShape([nodes(bcells(ic).int(ip).nen).x]',bcells(ic).int(ip).x,mCon.dm,mCon.pn);
            for neni=1:length(bcells(ic).int(ip).nen)
                T(2*neni-1:2*neni)=phi(neni).*bcells(ic).int(ip).bv;
                en(2*neni-1:2*neni)=[2*bcells(ic).int(ip).nen(neni)-1; 2*bcells(ic).int(ip).nen(neni)];
            end
            f(en)=f(en)+T*bcells(ic).int(ip).w*bcells(ic).J;
        end  
%     elseif bcells(ic).Bc==3
%     elseif bcells(ic).Bc==4
    end
end
%Assembly G matrix and q vector collocatd boundary particles
if mCon.BCuType~=1
    for i=1:bcCon.BCunx
        q(i)=-bcCon.BCux(2,i);
        [phi_i,dphidx_i,dphidy_i]=MLSShape([nodes(nodes(bcCon.BCux(1,i)).nen).x]',nodes(bcCon.BCux(1,i)).x,mCon.dm,mCon.pn);
        nj=0;
        for j=nodes(bcCon.BCux(1,i)).nen
            nj=nj+1;
            G(2*j-1,i)=G(2*j-1,i)-phi_i(nj);
        end
    end
    for i=1:bcCon.BCuny
        q([i+bcCon.BCunx])=-bcCon.BCuy(2,i);
        [phi_i,dphidx_i,dphidy_i]=MLSShape([nodes(nodes(bcCon.BCuy(1,i)).nen).x]',nodes(bcCon.BCuy(1,i)).x,mCon.dm,mCon.pn);
        nj=0;
        for j=nodes(bcCon.BCuy(1,i)).nen
            nj=nj+1;
            G(2*j,i+bcCon.BCunx)=G(2*j,i+bcCon.BCunx)-phi_i(nj);
        end
    end
end
time1=toc; %Assembly timer
tic %Solve timer
G = sparse(G);
end

