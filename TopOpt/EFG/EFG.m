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
% * _K_: the minimum stiffness matrix
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
    EFG(Ke,f,G,q,K,distrType,H,Hs,computeDerivatives)

GlobalConst
if nargin < 5
    [Ke,f,G,q,K] = EFGUnitMatrices();
end
if nargin < 6
    distrType = 1;
end
if nargin < 7 || isempty(H)
    filter = false;
else
    filter = true;
end
if nargin < 9
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
if computeDerivatives
    [dKdx{1:nd*mmCon.n, 1}] = deal(sparse(2*mCon.n,2*mCon.n));
end
dCdx = zeros(nd*mmCon.n,1);
mTot = 0;

if filter
    % If there is a filter, the densities are first computed and saved,
    % then filtered, then the matrices are assembled.
    massNodes = 1:nd*mmCon.n;
    rhoVec = zeros(mCon.m*mCon.nG^2,1);
    drhoVec = zeros(mCon.m*mCon.nG^2,nd*mmCon.n);
end


for ic=1:mCon.m                                         % Iterations over the internal cells
    for ip=1:cells(ic).ni                               % Iterations over the cell Gauss points
        if ~isempty(cells(ic).int(ip).nen) && ~isempty(cells(ic).int(ip).nemn)
            en=zeros(1,2*length(cells(ic).int(ip).nen));
            en(1:2:end-1)=2*[cells(ic).int(ip).nen]-1;      % x index of neighboring cells
            en(2:2:end)=2*[cells(ic).int(ip).nen];          % y index
            emn=zeros(1,nd*length(cells(ic).int(ip).nemn));
            for i = 1:nd
                emn(i:nd:end-nd+i)=nd*[cells(ic).int(ip).nemn]-nd+i;
            end
            [rho,drhodx] = asymptoticDensity(cells(ic).int(ip).x,...
                [mnodes(cells(ic).int(ip).nemn).x],...
                [mnodes(cells(ic).int(ip).nemn).theta],...
                [mnodes(cells(ic).int(ip).nemn).l]/2,...
                [mnodes(cells(ic).int(ip).nemn).m],...
                mmCon.rhoMax,distrType,true,mmCon.rm);
            if filter
                rhoVec((ic-1)*mCon.nG^2+ip) = rho;
            else
                K(en,en)=K(en,en)+(pCon.E-mmCon.EMin)*rho^oCon.p*Ke{(ic-1)*mCon.nG^2+ip};
            end
            mTot = mTot + rho*cells(ic).J*cells(ic).int(ip).w;
            if computeDerivatives
                if filter
                    drhoVec((ic-1)*mCon.nG^2+ip,emn) = drhodx;
                else
                    for i=1:length(emn)
                        dKdx{emn(i)}(en,en)=dKdx{emn(i)}(en,en)+...
                            (pCon.E-mmCon.EMin)*oCon.p*drhodx(i)*rho^(oCon.p-1)*Ke{(ic-1)*mCon.nG^2+ip};  
                    end
                end
            end
        end
    end
end

if filter
    % Filtering
    drhoVec = sparse(drhoVec);
    rhoVec = (H*rhoVec)./Hs;
    drhoVec = (H*drhoVec)./Hs;
    
    % Assembly
    for i = 1 : length(rhoVec)
        ic = ceil(i/mCon.nG^2);
        ip = mod(i-1,mCon.nG^2)+1;
        en=zeros(1,2*length(cells(ic).int(ip).nen));
        en(1:2:end-1)=2*[cells(ic).int(ip).nen]-1;      % x index of neighboring cells
        en(2:2:end)=2*[cells(ic).int(ip).nen];          % y index
        K(en,en)=K(en,en)+(pCon.E-mmCon.EMin)*rhoVec(i)^oCon.p*Ke{(ic-1)*mCon.nG^2+ip};
        if computeDerivatives
            emn = massNodes(drhoVec ~= 0);
            for j=1:length(emn)
                dKdx{emn(j)}(en,en)=dKdx{emn(j)}(en,en)+...
                    sparse((pCon.E-mmCon.EMin)*oCon.p*drhoVec(i,emn(j))*rhoVec(i)^(oCon.p-1)*Ke{(ic-1)*mCon.nG^2+ip});  
            end
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
if computeDerivatives
    for i = 1 : nd*mmCon.n
       dCdx(i) = -r(1:end-length(q))'*dKdx{i}*r(1:end-length(q));
    end
end
time3=toc;    

end
