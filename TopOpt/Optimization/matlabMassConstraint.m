% Matlab Mass constraint
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Similar to _massConstraint_ but adapted to the Matlab's algorithms

function [c,ceq,GC,GCeq] = matlabMassConstraint(x)
    
    [c,GC] = massConstraint(x);
    c = -c;
    GC = -GC;
    ceq = [];
    GCeq = [];

end