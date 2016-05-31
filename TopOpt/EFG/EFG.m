%% Element-Free Galerkin (EFG)
% Created by J.T.B. Overvelde - 18 April 2012
% 
% Master's thesis - The Moving Node Approach in Topology Optimization
% 
% <http://www.overvelde.com>
%
% Modified by G. Raze - June 2016
% 
% Determines the stiffness matrix for given nodal distribution and
% background mesh and determines the displacement by solving the system of
% equations for given boundary conditions.

function [ug,Compliance,dCdx,mTot,tk,tm,tdK,time1,time2,time3]=...
    EFG(Ke,f,G,q,distrType)

GlobalConst
if nargin < 4
    [Ke,f,G,q] = EFGUnitMatrices();
end
if nargin < 5
    distrType = 1;
end
tic %Assembly timer

%Assembly K matrix and f vector internal cells
cells = neighboringMassNodes(mnodes,cells);
if distrType == 3
    nd = 5;
else
    nd = distrType+1;
end
K=zeros(2*mCon.n);
[dKdx{1:nd*mmCon.n, 1}] = deal(sparse(2*mCon.n,2*mCon.n));
dCdx = zeros(nd*mmCon.n,1);
mTot = 0;
tk = 0;
tm = 0;
tdK = 0;
for ic=1:mCon.m                                         % Iterations over the internal cells
    for ip=1:cells(ic).ni                               % Iterations over the cell Guass points
        if ~isempty(cells(ic).int(ip).nen)
            en=zeros(1,2*length(cells(ic).int(ip).nen));
            en(1:2:end-1)=2*[cells(ic).int(ip).nen]-1;      % x index of neighboring cells
            en(2:2:end)=2*[cells(ic).int(ip).nen];          % y index
            emn=zeros(1,nd*length(cells(ic).int(ip).nemn));
            for i = 1:nd
                emn(i:nd:end-nd+i)=nd*[cells(ic).int(ip).nemn]-nd+i;
            end
            [rho,drhoAdx] = asymptoticDensity(cells(ic).int(ip).x,...
                [mnodes(cells(ic).int(ip).nemn).x],...
                [mnodes(cells(ic).int(ip).nemn).theta],...
                [mnodes(cells(ic).int(ip).nemn).l]/2,...
                [mnodes(cells(ic).int(ip).nemn).m],...
                mmCon.rhoMin,mmCon.rhoMax,distrType);
            tic
            K(en,en)=K(en,en)+rho*Ke{(ic-1)*mCon.nG^2+ip}(en,en);
            tk = tk+toc;
            tic
            mTot = mTot + rho*cells(ic).J*cells(ic).int(ip).w;
            tm = tm+toc;
            tic
            for i=1:length(emn)
                dKdx{emn(i)}(en,en)=dKdx{emn(i)}(en,en)+...
                    sparse(drhoAdx(i)*Ke{(ic-1)*mCon.nG^2+ip}(en,en)); 
            end
            tdK = tdK + toc;
        end
    end
end


time1=toc; %Assembly timer
%disp([num2str(time1),' seconds to assemble the matrices'])

tic %Solve timer

%Solve for the displacement constants belonging to the shape functions
K=sparse(K);

r=[K G;G' zeros(length(q))]\[f;q];
ug(:,1)=r(1:2:end-length(q)-1);
ug(:,2)=r(2:2:end-length(q));

%clear r K G f q

time2=toc; %Solve timer
%disp([num2str(time2),' seconds to solve the system'])
tic %Find nodal values timer

%Calculate the strains and stress constants belonging to the shape
%functions at the nodal positions
% eg=zeros(mCon.n,3); sg=zeros(mCon.n,3);
% for i=1:mCon.n
%     [phi,dphidx,dphidy]=MLSShape([nodes(nodes(i).nen).x]',nodes(i).x,mCon.dm,mCon.pn);
%     %strain
%     eg(i,1)=sum(dphidx.*ug([nodes(i).nen],1)'); 
%     eg(i,2)=sum(dphidy.*ug([nodes(i).nen],2)');
%     eg(i,3)=1/2*(sum(dphidx.*ug([nodes(i).nen],2)')...
%                  +sum(dphidy.*ug([nodes(i).nen],1)'));
%     %stress
%     sg(i,:)=pCon.D*eg(i,:)';
%     sg(i,3)=2*sg(i,3);
% end

time3=toc; %Find nodal values timer
%disp([num2str(time3),' seconds to find the nodal values'])
tic %Find error norm
 
Compliance=f'*r(1:end-length(q));
for i = 1 : nd*mmCon.n
   dCdx(i) = -r(1:end-length(q))'*dKdx{i}*r(1:end-length(q));
end
time4=toc;    
%disp([num2str(time4),' seconds to find the compliance'])


