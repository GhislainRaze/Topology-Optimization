%% Broyden-Fletcher-Goldfarb-Shanno quasi-Newton method
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Searches the minimum compliance with a quasi-Newton
% Broyden-Fletcher-Golfarb-Shanno algorithm. 
% The inputs are
%
% * _distrType_: the material distribution type (1, 2 or 3)
% * _method_: the discretization method type (1 or 2)

function history = quasiNewtonBFGS(distrType,method)

    preOptimization;
    I = eye(length(dCdx0));
    S0 = I;
    
    while relDif > oCon.relTol || deltaX > oCon.xTol && iter < oCon.iterMax
        iter = iter+1;

        x1 = x0 +s0*oCon.dg;
        if distrType >= 3
            x1 = checkFeasability(x1,x0);
        end
        if oCon.trueMinimum
            [x0p,C0p,dCdx0p,u0] = trueMinimum(x0,C0,dCdx0,x1,...
                objectiveFunction,oCon.iterMinimum, oCon.tolMinimum,...
                distrType>=3);
        else
            [x0p,C0p,dCdx0p,u0] = wolfe(x0,C0,dCdx0,x1,...
                objectiveFunction,oCon.iterWolfe,distrType>=3);
        end
        
        relDif = (C0-C0p)/C0p;

        y = (dCdx0p-dCdx0);
        gamma = 1/(s0'*y);
        S0 = (I-gamma*(s0*y'))*S0*(I-gamma*(y*s0')) + gamma*s0*s0';
        s0 = -S0*dCdx0p;
        s0 = s0/norm(s0);
        postIteration;
    end
    
    postOptimization;
end