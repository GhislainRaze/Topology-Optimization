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


function [x1,x1Changed] = checkFeasability(x1)

    global oCon mmCon
    
    x1Changed = false;
    [cm,dcmdx] = massConstraint(x1,oCon.relaxation);
    if cm < 0
        x1Changed = true;
        nm = length(x1)/5;
        ind = zeros(2*nm,1);
        for i = 1 : nm
            ind(2*i-1) = 5*i-1;
            ind(2*i) = 5*i;
        end
        H = mmCon.rm/4*sparse(kron(eye(nm),[0 1; 1 0]));
        
        a = 0.5*dcmdx(ind)'*H*dcmdx(ind);
        b = dcmdx(ind)'*H*x1(ind);
        
        alpha = (-b - sqrt(b^2+4*a*cm))/(2*a);
        
        x1(ind) = x1(ind) + 1.1*alpha*dcmdx(ind);
    end
end