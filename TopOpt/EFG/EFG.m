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
% * _H_: the filter convolution matrix (optional)
% * _Hs_: the sum of the filter convolution matrix lines (optional)
% *_computeDerivatives_: false if the compliance derivatives do not need to
% be computed (default value = true)
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
% * _time1_: the time to make the assembly
% * _time2_: the time to solve the linear system
% * _time3_: the time to compute the compliance and its sensitivities


function [ug,Compliance,dCdx,mTot,time1,time2,time3]=...
    EFG(Ke,f,G,q,distrType,H,Hs,computeDerivatives)

GlobalConst
if nargin < 4
    [Ke,f,G,q] = EFGUnitMatrices();
end
if nargin < 5
    distrType = 1;
end
if nargin < 6 || isempty(H)
    filter = false;
else
    filter = true;
end
if nargin < 8
    computeDerivatives = true;
end
tic %Assembly timer

%Assembly K matrix and f vector internal cells
cells = neighboringMassNodes(mnodes,cells);
if distrType == 3
    nd = 5;
else
    nd = distrType+1;
end

rhoVec = zeros(mCon.m*mCon.nG^2,1);
drhoVec = zeros(mCon.m*mCon.nG^2,nd*mmCon.n);
K = zeros(2*mCon.n);

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
            [rhoVec((ic-1)*mCon.nG^2+ip),drhoVec((ic-1)*mCon.nG^2+ip,emn)] = asymptoticDensity(cells(ic).int(ip).x,...
                mnodes(cells(ic).int(ip).nemn),pCon.filledRegions,mmCon.rf,...
                mmCon.rm,mmCon.rhoMax,distrType,true);
            if ~filter
                K(en,en) = K(en,en) + (mmCon.EMin+(pCon.E-mmCon.EMin)*rhoVec((ic-1)*mCon.nG^2+ip)^oCon.p)*Ke{(ic-1)*mCon.nG^2+ip};
            end
        end
    end
end
mTot = rhoVec'*mCon.w;

if filter
    % Filtering
    rhoVec = (H*rhoVec)./Hs;
    for i = 1 : size(drhoVec,2)
        drhoVec(:,i) = H*(drhoVec(:,i)./Hs);
    end
    EVec = mmCon.EMin+(pCon.E-mmCon.EMin).*rhoVec.^oCon.p;
    for i = 1 : length(rhoVec)
        ic = ceil(i/mCon.nG^2);
        ip = mod(i-1,mCon.nG^2)+1;
        en=zeros(1,2*length(cells(ic).int(ip).nen));
        en(1:2:end-1)=2*[cells(ic).int(ip).nen]-1;      % x index of neighboring cells
        en(2:2:end)=2*[cells(ic).int(ip).nen];          % y index
        K(en,en)=K(en,en)+EVec(i)*Ke{i};
    end
end
drhoVec = sparse(drhoVec);
K = (K+K')/2;
K = sparse(K);

time1=toc; %Assembly timer

tic %Solve timer

%Solve for the displacement constants belonging to the shape functions

r = [K G;G' zeros(length(q))]\[f;q];
ug(:,1)=r(1:2:end-length(q)-1);
ug(:,2)=r(2:2:end-length(q));


time2=toc; %Solve timer

tic %Find error norm

if pCon.type == 1
    Compliance=f'*r(1:end-length(q));
elseif pCon.type == 2
    Compliance=f(:,2)'*r(1:end-length(q),1);
end
if computeDerivatives
    Ce = zeros(mCon.m*mCon.nG^2,1);
    for ic = 1 : mCon.m
        for ip=1:cells(ic).ni
            en=zeros(1,2*length(cells(ic).int(ip).nen));
            en(1:2:end-1)=2*[cells(ic).int(ip).nen]-1;
            en(2:2:end)=2*[cells(ic).int(ip).nen];
            if pCon.type == 1
                Ce((ic-1)*mCon.nG^2+ip) = r(en)'*Ke{(ic-1)*mCon.nG^2+ip}*r(en);
            elseif pCon.type == 2
                Ce((ic-1)*mCon.nG^2+ip) = r(en,2)'*Ke{(ic-1)*mCon.nG^2+ip}*r(en,1);
            end
        end
    end
    dCdx = -(diag((pCon.E-mmCon.EMin)*oCon.p*rhoVec.^(oCon.p-1))*drhoVec)'*Ce;
else
    dCdx = zeros(nd*mmCon.n,1);   
end
time3=toc;    

end
