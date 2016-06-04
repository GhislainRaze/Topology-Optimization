%% Steepest descent method
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Initialize problem and mesh for a FEM analysis. The global constants are
% modified.

function steepestDescent(distrType,method)
    GlobalConst
    iter = 0;
    relDif = 1;
    if method == 1
        [Ke,f,G,q]=EFGUnitMatrices();
        objfun = @(x) complianceEFG(x,distrType,Ke,f,G,q);
    elseif method == 2
        [Ke,f,ubar]=FEMUnitMatrices();
        objfun = @(x) complianceFEM(x,distrType,Ke,f,ubar);
    end
    x0 = mnodesToVector(mnodes,distrType);
    [C0,dCdx0] = objfun(x0);
    s0 = - dCdx0/norm(dCdx0);

    while abs(relDif) > oCon.relTol && iter < oCon.iterMax
        iter = iter+1;

        x1 = x0 - dCdx0/norm(dCdx0)*oCon.dg;

        if oCon.trueMinimum
            [x0,C0p,dCdx0] = trueMinimum(x0,C0,dCdx0,x1,objfun);
        else
            [x0,C0p,dCdx0] = wolfe(x0,C0,dCdx0,x1,objfun);
        end

        relDif = (C0-C0p)/C0p;
        C0 = C0p;

        disp(['Iteration ',num2str(iter),' : Compliance = ',...
            num2str(C0), ' ; Gradient norm = ',num2str(norm(dCdx0))])
    end

end
