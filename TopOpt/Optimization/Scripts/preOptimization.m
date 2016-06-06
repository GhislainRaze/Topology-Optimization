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

% Objective function definition
if method == 1
    [Ke,f,G,q]=EFGUnitMatrices();
    if oCon.filter && ~oCon.filterIter 
        [H,Hs] = filterInitialization(cells,mCon.nG,mmCon.rmin);
        objectiveFunction = @(x) complianceEFG(x,distrType,Ke,f,G,q,H,Hs);
        filterEnabled = true;
        disp('Filter enabled')
    else
        objectiveFunction = @(x) complianceEFG(x,distrType,Ke,f,G,q);
        filterEnabled = false;
    end
elseif method == 2
    [Ke,f,ubar]=FEMUnitMatrices();
    if oCon.filter && ~oCon.filterIter 
        [H,Hs] = filterInitialization(cells,mCon.nG,mmCon.rmin);
        objectiveFunction = @(x) complianceFEM(x,distrType,Ke,f,ubar,H,Hs);
        filterEnabled = true;
        disp('Filter enabled')
    else
        objectiveFunction = @(x) complianceFEM(x,distrType,Ke,f,ubar);
        filterEnabled = false;
    end
end

% Variables initialization
x0 = mnodesToVector(mnodes,distrType);
[C0,dCdx0] = objectiveFunction(x0);
s0 = - dCdx0/norm(dCdx0);

% Display information
disp(['Iteration ',num2str(iter),' : Compliance = ',...
    num2str(C0), ' ; Gradient norm = ',num2str(norm(dCdx0))])

% History initialization
history.x = x0;
history.C = C0;