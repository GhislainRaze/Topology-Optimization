%% Pre Optimization steps
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Initializes the optimization variables.

GlobalConst

% Initializations
iter = 0;
relDif = 1;
deltaX = 1;
time = 0;


% Objective function definition
if method == 1
        [Ke,f,G,q,K]=EFGUnitMatrices();
        disp('Unit matrices computed')
        if oCon.filter && ~oCon.filterIter 
            [H,Hs] = filterInitialization(cells,mCon.nG,mmCon.rmin);
            filterEnabled = true;
            disp('Filter enabled')
        else
            H = [];
            Hs = [];
            filterEnabled = false;
        end
        objectiveFunction = @(x) complianceEFG(x,distrType,Ke,f,G,q,K,...
                H,Hs);
    elseif method == 2
        [Ke,f,ubar,K]=FEMUnitMatrices();
        disp('Unit matrices computed')
        if oCon.filter && ~oCon.filterIter 
            [H,Hs] = filterInitialization(cells,mCon.nG,mmCon.rmin);
            filterEnabled = true;
            disp('Filter enabled')
        else
            H = [];
            Hs = [];
            filterEnabled = false;
        end
        objectiveFunction = @(x) complianceFEM(x,distrType,Ke,f,ubar,K,...
                H,Hs);
    end

% Check mesh and mass distribution
volFrac = mmCon.vol/pCon.vol;
disp(['Volume fraction: ',num2str(100*volFrac),'%'])
a = min(mmCon.dx,mmCon.dy)/max(mCon.dx,mCon.dy);
if a < 1
    disp('Warning: the mesh is probably too coarse for the mass distribution.')
    disp(['    The ratio between mass nodes influence domain and cells/element',...
        ' size should be at least 1'])
    disp(['    Current value: ',num2str(a)])
end

% Variables initialization
x0 = mnodesToVector(mnodes,distrType);
[C0,dCdx0,u0] = objectiveFunction(x0);
s0 = - dCdx0/norm(dCdx0);


% Display information
disp(['Iteration ',num2str(iter),' : Compliance = ',...
    num2str(C0), ' ; Gradient norm = ',num2str(norm(dCdx0))])

% History initialization
history.x = x0;
history.C = C0;
history.u = u0;
history.m = mTot;

tic2 = tic;