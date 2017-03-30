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


function history = BFGSMax(distrType,method)

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
    nSuppressed = 1;
    
    while nMerged > 0 || nSuppressed > 0;
        
        
        I = eye(length(dCdx0));
        S0 = I;
        while abs(relDif) > oCon.relTol && deltaX > oCon.xTol && iter < oCon.iterMax
            iter = iter+1;
            
            if homogenous
            for i = 1:nd:length(dCdx0)
                dCdx0(i:i+nd-1) = dt*dCdx0(i:i+nd-1)/norm(dCdx0(i:i+nd-1));
            end
            end
            
            s0 = -S0*dCdx0;
            x0p = x0 + dt*s0/norm(s0);

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
            
            if length(x0) ~= length(x0p)
                disp(num2str(length(x0)))
                disp(num2str(length(x0p)))
            end
            
            if length(dCdx0) ~= length(dCdx0p)
                disp(num2str(length(dCdx0)))
                disp(num2str(length(dCdx0p)))
            end
            
            s0 = (x0p-x0)/norm(x0p-x0);
            y = dCdx0p-dCdx0;
            gamma = 1/(s0'*y);
            S0 = (I-gamma*(s0*y'))*S0*(I-gamma*(y*s0')) + gamma*(s0*s0');

            postIteration;



        end

        postConvergence
        
    end
    
    postOptimization;
end