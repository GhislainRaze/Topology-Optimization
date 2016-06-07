%% True minimum linesearch
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Searches the true minimum of a function along direction _x1-x0_, starting
% from _x0_. The minimum is considered as a true minimum if $f'(x^*) <
% tol\times f'(x_0)$.
%
% The inputs are
%
% * _x0_: the starting point
% * _f0_: the function value at _x0_
% * _g0_: the gradient value at _x0_
% * _x1_: a point indicating the search direction (_x1-x0_)
% * _objfun_: the objective function
% * _iterMax_: the maximum number of iterations (optional, default value =
% 100)
% * _tol_: the tolerance on the minimum derivative (optional, default value
% = 0.001)

function  [xs,fs,gs,us] = trueMinimum(x0,f0,g0,x1,objfun,iterMax,tol)

    if nargin < 6
        iterMax = 100;
    end
    if nargin < 7
        tol = 1e-3;
    end
    minimum = false;
    iter = 0;
    s0 = (x1-x0)/norm(x1-x0);
    
    [f1,g1] = objfun(x1);
    
    gp0 = g0'*s0;
    gp1 = g1'*s0;
    g00 = norm(g0);
    
    while ~minimum && iter < iterMax
        
        iter = iter+1;
        
        ls = cubicApproximation(0,norm(x1-x0),f0,f1,gp0,gp1);
        
        xs = x0+ls*s0;
        [fs,gs,us] = objfun(xs);
        gps = gs'*s0;
        
        if abs(gps/g00) > tol
            if norm(x1-xs) > norm(xs-x0)
                x0 = xs;
                f0 = fs;
                gp0 = gps;
            else
                x1 = xs;
                f1 = fs;
                gp1 = gps;
            end
        else
            minimum = true;
        end
    end
end