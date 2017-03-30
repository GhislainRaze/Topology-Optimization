%% Overvelde's Optimization Algorithm
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Calls Overvelde's algorithm to optimize the material distribution.
%
% Overvelde's algorithm is a flow-like optimization algorithm based on a
% steepest descend and a fixed step.
%
% The flow-like equation
%
% $$ \frac{\partial^2 x_i^I(t)}{\partial t^2} + c \frac{\partial
% x_i^I(t)}{\partial t} = -\frac{\partial C(x_i^I(t))}{\partial x_i^I(t)} $$
%
% is discretized in the following way (Euler explicit time integration)
%
% $$ v_i^I(t^{k+1}) = v_i^I(t^k) - \Delta t \left(c v_i^I(t^k) + \frac{\partial 
% C(x_i^I(t^k))}{\partial x_i^I(t^k)}\right)$$
%
% $$ x_i^I(t^{k+1}) = x_i^I(t^k) + \Delta t v_i^I(t^{k+1})$$


function history = gradientMax(distrType,method)

    preOptimization;
    moveLimit = oCon.dg;
    dt = oCon.dg;
    dtmax = dt;
    dtmin = 0;
    homogenous = true;
    
    m1 = 0.01;
    m2 = 0.15;
    
    
    if distrType == 3
        nd = 5;
    else
        nd = distrType+1;
    end
    
    nMerged = 1;
    
    while nMerged > 0
        
        
        while abs(relDif) > oCon.relTol && deltaX > oCon.xTol && iter < oCon.iterMax
            iter = iter+1;
            if homogenous
            for i = 1:nd:length(dCdx0)
                dCdx0(i:i+nd-1) = dt*dCdx0(i:i+nd-1)/norm(dCdx0(i:i+nd-1));
            end
            else
                dCdx0 = dt*dCdx0/norm(dCdx0);
            end
            x0p = x0 - dCdx0;

            if distrType >= 3
                x0p = checkFeasability(x0p);
            end
            [C0p,dCdx0p,u0] = objectiveFunction(x0p);

            if C0p > C0+m1*dt*norm(dCdx0)      % First Wolfe criterion
                dt = (dtmin+dt)/2;
                dtmax = dt;
                dtmin = dt/2;
            elseif norm(dCdx0p) < m2*norm(dCdx0)        % Second Wolfe crietrion
                dt = (dt+dtmax)/2;
                dtmin = dt;
                dtmax = min(2*dt,moveLimit);
            end


            postIteration;



        end


        postConvergence
    
    end
    
    postOptimization;
end