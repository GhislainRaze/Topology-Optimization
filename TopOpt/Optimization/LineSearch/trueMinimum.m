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
% * _objectiveFunction_: the objective function
% * _iterMax_: the maximum number of iterations (optional, default value =
% 100)
% * _tol_: the tolerance on the minimum derivative (optional, default value
% = 0.001)

function  [xs,fs,gs,us] = trueMinimum(x0,f0,g0,x1,objectiveFunction,...
    iterMax,tol,massC)

    if nargin < 6
        iterMax = 100;
    end
    if nargin < 7
        tol = 1e-3;
    end
    minimum = false;
    iter = 0;
    s0 = (x1-x0)/norm(x1-x0);
    
    [f1,g1] = objectiveFunction(x1);
    
    gp0 = g0'*s0;
    gp1 = g1'*s0;
    g00 = norm(g0);
    
    while ~minimum && iter < iterMax
        
        iter = iter+1;
        
        ls = cubicApproximation(0,norm(x1-x0),f0,f1,gp0,gp1);
        
        xs = x0+ls*s0;
        if massC
            [xs,changedDirection] = checkFeasability(xs,x0);
            [fs,gs,us] = objectiveFunction(xs);
            if changedDirection
                s0 = (xs-x0)/norm(xs-x0);
                gp0 = g0'*s0;
                gp1 = g1'*s0;
            end
        else
            [fs,gs,us] = objectiveFunction(xs);
        end
        
        gps = gs'*s0;
        
        if abs(gps/g00) > tol
            if gps > 0
               x1 = xs;
               f1 = fs;
               gp1 = gps;
               x0 = (x0+x1)/2;
               [f0,g0] = objectiveFunction(x0);
               gp0 = g0'*s0;
            else
               x0 = xs;
               f0 = fs;
               gp0 = gps;
               x1 = (x0+x1)/2;
               [f1,g1] = objectiveFunction(x1);
               gp1 = g1'*s0;
            end
        else
            minimum = true;
        end
    end
end