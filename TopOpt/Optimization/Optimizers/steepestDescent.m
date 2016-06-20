%% Steepest descent method
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Searches the minimum compliance with a steepest descent algorithm.
% The inputs are
%
% * _distrType_: the material distribution type (1, 2 or 3)
% * _method_: the discretization method type (1 or 2)

function history = steepestDescent(distrType,method)
    preOptimization;
    step = oCon.dg;
    
    while relDif > oCon.relTol && deltaX > oCon.xTol && iter < oCon.iterMax
        iter = iter+1;

        x1 = x0 - dCdx0*step;

        if distrType >= 3
        	x1 = checkFeasability(x1);
        end
        if oCon.trueMinimum
            [x0p,C0p,dCdx0p,u0] = trueMinimum(x0,C0,dCdx0,x1,...
                objectiveFunction,oCon.iterMinimum, oCon.tolMinimum,...
                distrType>=3);
        else
            [x0p,C0p,dCdx0p,u0] = wolfe(x0,C0,dCdx0,x1,...
                objectiveFunction,oCon.iterWolfe,distrType>=3);
        end
        
        step = min(1,2.02*(C0p-C0)/norm(dXdc0p));
        postIteration;
    end
postOptimization;
end
