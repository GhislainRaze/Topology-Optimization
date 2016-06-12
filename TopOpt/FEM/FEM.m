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
% * _K_: the minimum stiffness matrix
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
    FEM(Ke,f,ubar,K,distrType,H,Hs,computeDerivatives)

GlobalConst
if nargin < 4
    [Ke,f,ubar,K] = FEMUnitMatrices();
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

%% Assembly K matrix and f vector internal cells
tic %Assembly timer
for ic=1:mCon.m                                         % Iterations over the internal cells
    en=zeros(1,2*length(cells(ic).nen));
    en(1:2:end-1)=2*[cells(ic).nen]-1;              % x index of neighboring cells
    en(2:2:end)=2*[cells(ic).nen];                  % y index
    for ip=1:cells(ic).ni                               % Iterations over the cell Gauss points
        if ~isempty(cells(ic).int(ip).nemn)
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
    for i = 1 : size(drhoVec,2)
        drhoVec(:,i) = H*(drhoVec(:,i)./Hs);
    end
    
    % Assembly
    for i = 1 : length(rhoVec)
        if ~rhoVec(i)
            ic = ceil(i/mCon.nG^2);
            ip = mod(i-1,mCon.nG^2)+1;
            en=zeros(1,2*length(cells(ic).nen));
            en(1:2:end-1)=2*[cells(ic).nen]-1;      % x index of neighboring cells
            en(2:2:end)=2*[cells(ic).nen];          % y index
            K(en,en)=K(en,en)+(pCon.E-mmCon.EMin)*rhoVec(i)^oCon.p*Ke{(ic-1)*mCon.nG^2+ip};
            if computeDerivatives
                emn=zeros(1,nd*length(cells(ic).int(ip).nemn));
                for j = 1:nd
                    emn(j:nd:end-nd+j)=nd*[cells(ic).int(ip).nemn]-nd+j;
                end
                for j=1:length(emn)
                    dKdx{emn(j)}(en,en)=dKdx{emn(j)}(en,en)+...
                        (pCon.E-mmCon.EMin)*oCon.p*drhoVec(i,emn(j))*rhoVec(i)^(oCon.p-1)*Ke{(ic-1)*mCon.nG^2+ip}; 
                end
            end
        end
   end
end

time1=toc; %Assembly timer


%% Solver
tic %Solve timer
% Partitioning between the known and unknown displacements
cn = find(~isnan(ubar));        % Constrained nodes indices
fn = 1:2*mCon.n;                % Free nodes indices
fn(cn) =[];
Kic = K(fn,fn);
Kcc = K(fn,cn);
uc = ubar(cn);
fc = f(fn);



u=Kic\(fc-Kcc*uc);             % Solution of the linear system for the free nodes
r = zeros(2*mCon.n,1);
r(cn) = ubar(cn);               % Imposed nodal displacements
r(fn) = u;                      % Computed nodal displacements

ug(:,1)=r(1:2:end-1);
ug(:,2)=r(2:2:end);



time2=toc; %Solve timer

 
%% Computing the compliance and its derivatives

tic 
Compliance=f'*r;
if computeDerivatives
    for i = 1 : nd*mmCon.n
       dCdx(i) = -r'*dKdx{i}*r;
    end
end
time3=toc;    
end

