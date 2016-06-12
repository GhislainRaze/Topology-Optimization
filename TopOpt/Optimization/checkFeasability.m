%% Check feasability
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Checks if the mass nodes variables vector _x1_ represents a feasible
% configuration (ie if the structural mass does not exceed the maximum
% allowed mass). If not, the constraint gradient is subtracked to the
% vector _x1-x0_ up until _x1_ represents a feasible configuration.
%
% This function returns a feasible _x1_ and a flag _x1Changed_ which is
% true if _x1_ has changed.


function [x1,x1Changed] = checkFeasability(x1,x0)

    global oCon
    
    x1Changed = false;
    
    % Compute the mass constraint and gradient
    [cm,dcmdx] = massConstraint(x1,oCon.relaxation);
    
    % Subtrack the constraint gradient
    if cm <=0
        step0 = 1.01*dcmdx'*(x1-x0)/(norm(dcmdx)^2);
        step = step0;
        x1Changed = true;
        while cm <= 0
            x1 = x1 - step*dcmdx;
            if step > 1e-10
            	cm = massConstraint(x1,oCon.relaxation);
            else
                [cm,dcmdx] = massConstraint(x1,oCon.relaxation);
            end
            step = 0.01*step0;
        end
    end
end