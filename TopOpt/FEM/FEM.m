%% Finite Element Method (FEM)
%
%

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
        en(1:2:end-1)=2*[cells(ic).nen]-1;      % x index of neighboring cells
        en(2:2:end)=2*[cells(ic).nen];          % y index
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


time1=toc; %Assembly timer
%disp([num2str(time1),' seconds to assemble the matrices'])


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


%clear r K G f q

time2=toc; %Solve timer
%disp([num2str(time2),' seconds to solve the system'])
tic %Find nodal values timer

time3=toc; %Find nodal values timer
%disp([num2str(time3),' seconds to find the nodal values'])

 
%% Computing the compliance and its derivatives

tic 
Compliance=f'*r;
for i = 1 : nd*mmCon.n
   dCdx(i) = -r'*dKdx{i}*r;
end
time4=toc;    
%disp([num2str(time4),' seconds to find the compliance'])
end

