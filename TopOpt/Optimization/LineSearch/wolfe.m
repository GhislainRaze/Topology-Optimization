%% Wolfe minimum linesearch
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Searches the minimum of a function. The minimum is not exactly located,
% but is deemed satisfactory if the two Wolfe criteria are satisfied.
%
% The first Wolfe criterion states that the function value must be low
% enough compared to that of the starting point:
%
% $$ \phi(x) \leq \phi (0) + m_1 \phi'(0)x $$
%
% with $m_1 = 0.01$. The second Wolfe criterion states that the function
% derivative must be increased enough compared to that of the starting 
% point:
%
% $$ \phi'(x) \geq m_2 \phi'(0) $$
%
% with $m_2 = 0.15$. The inputs are
%
% * _x0_ and _x1_: two starting points
% * _f0_: the function value at point _x0_
% * _g0_: the gradient at point _x0_


function  [xs,fs,gs] = wolfe(x0,f0,g0,x1,objfun)

    m1 = 0.01;
    m2 = 0.15;
    iterMax = 10;
    minimum = false;
    iter = 0;
    s0 = (x1-x0)/norm(x1-x0);
    
    [f1,g1] = objfun(x1);
    gp0 = g0'*s0;
    gp1 = g1'*s0;
    
    f00 = f0;
    gp00 = gp0;
    
    while ~minimum && iter < iterMax
        
        iter = iter + 1;
        
        ls = cubicApproximation(0,norm(x1-x0),f0,f1,gp0,gp1);
        xs = x0+ls*s0;
        
        [fs,gs] = objfun(xs);
        gps = gs'*s0;
        
        if fs > f00+m1*ls*gp00      % First Wolfe criterion
            x1 = xs;
            f1 = fs;
            g1 = gs;
            gp1 = g1'*s0;
        elseif gps < m2*gp00        % Second Wolfe crietrion
            x0 = xs;
            f0 = fs;
            g0 = gs;
            gp0 = g0'*s0;
        else
            minimum = true;
        end
    end
    
end