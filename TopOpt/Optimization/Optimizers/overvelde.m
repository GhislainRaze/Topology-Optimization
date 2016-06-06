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


function history = overvelde(distrType,method)

    preOptimization;
    c = 0.8;
    x0p = x0;
    dt = 1/(max(mCon.nx,mCon.ny));
    v = zeros(size(x0));
    while abs(relDif) > oCon.relTol && deltaX > oCon.xTol && iter < oCon.iterMax
            iter = iter+1;
        
        v = v - dt*(c*v+dCdx0/norm(dCdx0));
        x0p = x0 + v*dt;
        

        [C0p,dCdx0p] = objectiveFunction(x0p);
        postIteration;
    end

    postOptimization;
end