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
        
        a = dcmdx(ind)'*H*dcmdx(ind);
        b = 0.5*(dcmdx(ind)'*H*x1(ind)+x1(ind)'*H*dcmdx(ind));
        
        alpha1 = (-b + sqrt(b^2+2*a*cm))/a;
        alpha2 = (-b - sqrt(b^2+2*a*cm))/a;
        
        x11 = x1;
        x12 = x1;
        x11(ind) = x1(ind) + alpha1*dcmdx(ind);
        x12(ind) = x1(ind) + alpha2*dcmdx(ind);
        if min(x12(ind)) < 0
            x1 = x11;
        else
            x1 = x12;
        end
    end
end