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
% * _tK_: the time to assemble the stiffness matrix
% * _tm_: the time to compute the total mass
% * _tdK_: the time to assemble the derivative of the stiffness matrix
% * _time1_: the time to make the assembly
% * _time2_: the time to solve the linear system
% * _time3_: the time to compute the compliance and its sensitivities


function [ug,Compliance,dCdx,mTot,tk,tm,tdK,time1,time2,time3]=...
    FEM(Ke,f,ubar,distrType)

GlobalConst
if nargin < 3
    [Ke,f,ubar] = FEMUnitMatrices();
end
if nargin < 4
    distrType = 1;
end



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

%% Assembly K matrix and f vector internal cells
tic %Assembly timer
for ic=1:mCon.m                                         % Iterations over the internal cells
    for ip=1:cells(ic).ni                               % Iterations over the cell Gauss points
        en=zeros(1,2*length(cells(ic).nen));
        en(1:2:end-1)=2*[cells(ic).nen]-1;              % x index of neighboring cells
        en(2:2:end)=2*[cells(ic).nen];                  % y index
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


time1=toc; %Assembly timer


%% Solver
tic %Solve timer
K=sparse(K);
% Partitioning between the known and unknown displacements

cn = find(~isnan(ubar));        % Constrained nodes indices
fn = 1:2*mCon.n;                % Free nodes indices
fn(cn) =[];
Kic = K(fn,fn);
Kcc = K(fn,cn);
uc = ubar(cn);
fc = f(fn);

u=Kic\(fc-Kcc*uc);              % Solution of the linear system for the free nodes
r = zeros(2*mCon.n,1);
r(cn) = ubar(cn);               % Imposed nodal displacements
r(fn) = u;                      % Computed nodal displacements

ug(:,1)=r(1:2:end-1);
ug(:,2)=r(2:2:end);



time2=toc; %Solve timer

 
%% Computing the compliance and its derivatives

tic 
Compliance=f'*r;
for i = 1 : nd*mmCon.n
   dCdx(i) = -r'*dKdx{i}*r;
end
time3=toc;    
end

