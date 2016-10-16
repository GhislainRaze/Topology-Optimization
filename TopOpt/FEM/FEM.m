%% Finite Element Method (FEM)
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
% * _Ke_: the unit stiffness matrix
% * _f_: the nodal force vector
% * _ubar_: the imposed nodal displacements (set to NaN if there is no
% imposed displacement associated to a given node)
% * _distrType_: the mass distribution type
% * _H_: the filter convolution matrix (optional)
% * _Hs_: the sum of the filter convolution matrix lines (optional)
% *_computeDerivatives_: false if the compliance derivatives do not need to
% be computed (default value = true)
%
% If they are not specified, they are computed from FEMUnitMatrices
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
    FEM(Ke,f,ubar,distrType,H,Hs,computeDerivatives)

GlobalConst
if nargin < 3
    [Ke,f,ubar] = FEMUnitMatrices();
end
if nargin < 4
    distrType = 1;
end
if nargin < 5 || isempty(H)
    filter = false;
else
    filter = true;
end
if nargin < 7
    computeDerivatives = true;
end

cells = neighboringMassNodes(mnodes,cells);
if distrType == 3
    nd = 5;
else
    nd = distrType+1;
end


rhoVec = zeros(mCon.m*mCon.nG^2,1);
drhoVec = zeros(mCon.m*mCon.nG^2,nd*mmCon.n);


%% Assembly K matrix and f vector internal cells
tic %Assembly timer
for ic=1:mCon.m                                         % Iterations over the internal cells
    for ip=1:cells(ic).ni                               % Iterations over the cell Gauss points
        emn=zeros(1,nd*length(cells(ic).int(ip).nemn));
        for i = 1:nd
            emn(i:nd:end-nd+i)=nd*[cells(ic).int(ip).nemn]-nd+i;
        end
        [rhoVec((ic-1)*mCon.nG^2+ip),drhoVec((ic-1)*mCon.nG^2+ip,emn)] = asymptoticDensity(cells(ic).int(ip).x,...
            mnodes(cells(ic).int(ip).nemn),pCon.filledRegions,mmCon.rf,...
            mmCon.rm,mmCon.rhoMax,distrType,true);
    end
end

mTot = rhoVec'*mCon.w;

if filter
    % Filtering
    rhoVec = (H*rhoVec)./Hs;
    for i = 1 : size(drhoVec,2)
        drhoVec(:,i) = H*(drhoVec(:,i)./Hs);
    end
end
drhoVec = sparse(drhoVec);


nEDofs = (2*length(cells(1).nen))^2;

EVec = (mmCon.EMin+(pCon.E-mmCon.EMin).*rhoVec.^oCon.p);
sK = reshape(reshape(cell2mat(Ke(:))',nEDofs,mCon.nG^2)*...
    reshape(EVec,mCon.nG^2,mCon.m),mCon.m*nEDofs,1);
K = sparse(mCon.iK,mCon.jK,sK);
K = (K+K')/2;



time1=toc; %Assembly timer


%% Solver
tic %Solve timer
% Partitioning between the known and unknown displacements
cn = find(~isnan(ubar));        % Constrained nodes indices
fn = 1:2*mCon.n;                % Free nodes indices
fn(cn) =[];
r = zeros(2*mCon.n,size(f,2));


r(fn,:)=K(fn,fn)\(f(fn,:)-kron(K(fn,cn)*ubar(cn),ones(1,size(f,2))));             % Solution of the linear system for the free nodes
r(cn) = ubar(cn);               % Imposed nodal displacements

ug(:,1)=r(1:2:end-1);
ug(:,2)=r(2:2:end);



time2=toc; %Solve timer

 
%% Computing the compliance and its derivatives

tic 
if pCon.type == 1
    Compliance=f'*r;
elseif pCon.type == 2
    Compliance=f(:,2)'*r(:,1);
end
if computeDerivatives
    Ce = zeros(mCon.m*mCon.nG^2,1);
    for ic = 1 : mCon.m
        en=zeros(1,2*length(cells(ic).nen));
        en(1:2:end-1)=2*[cells(ic).nen]-1;              % x index of neighboring cells
        en(2:2:end)=2*[cells(ic).nen];                  % y index
        for ip=1:cells(ic).ni
            if pCon.type == 1
                Ce((ic-1)*mCon.nG^2+ip) = r(en)'*Ke{ip}*r(en);
            elseif pCon.type == 2
                Ce((ic-1)*mCon.nG^2+ip) = r(en,2)'*Ke{ip}*r(en,1);
            end
        end
    end
    dCdx = -(diag((pCon.E-mmCon.EMin).*rhoVec.^(oCon.p-1)*oCon.p)*drhoVec)'*Ce;
else
    dCdx = zeros(nd*mmCon.n,1);   
end
time3=toc;    
end

