%% Element-Free Galerkin (EFG)
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Determines the stiffness matrix for given nodal distribution, material 
% distribution and background mesh and determines the displacement by 
% solving the system of equations for given boundary conditions.
%
%
% The inputs are
%
% * _Ke_: the unit stiffness matrices
% * _f_: the nodal force vector
% * _G_: the lagrangian multipliers matrix
% * _q_: the second member associated with lagrangian multipliers
% * _distrType_: the material distribution type
%
% If they are not specified, they are computed from EFGUnitMatrices
%
% The outputs are
%
% * _ug_: the nodal displacements
% * _Compliance_
% * _dCdx_: the compliance sensitivities with respect to the material
% distribution variables
% * _mTot_: the total mass of the structure
% * _tK_: the time to assemble the stiffness matrix
% * _tm_: the time to compute the total mass
% * _tdK_: the time to assemble the derivative of the stiffness matrix
% * _time1_: the time to make the assembly
% * _time2_: the time to solve the linear system
% * _time3_: the time to compute the compliance and its sensitivities


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
    for ip=1:cells(ic).ni                               % Iterations over the cell Gauss points
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
            K(en,en)=K(en,en)+rho^mmCon.p*Ke{(ic-1)*mCon.nG^2+ip};
            tk = tk+toc;
            tic
            mTot = mTot + rho*cells(ic).J*cells(ic).int(ip).w;
            tm = tm+toc;
            tic
            for i=1:length(emn)
                dKdx{emn(i)}(en,en)=dKdx{emn(i)}(en,en)+...
                    sparse(mmCon.p*drhoAdx(i)*rho^(mmCon.p-1)*Ke{(ic-1)*mCon.nG^2+ip}); 
            end
            tdK = tdK + toc;
        end
    end
end


time1=toc; %Assembly timer

tic %Solve timer

%Solve for the displacement constants belonging to the shape functions
K=sparse(K);

r=[K G;G' zeros(length(q))]\[f;q];
ug(:,1)=r(1:2:end-length(q)-1);
ug(:,2)=r(2:2:end-length(q));


time2=toc; %Solve timer

tic %Find error norm
 
Compliance=f'*r(1:end-length(q));
for i = 1 : nd*mmCon.n
   dCdx(i) = -r(1:end-length(q))'*dKdx{i}*r(1:end-length(q));
end
time3=toc;    

end
