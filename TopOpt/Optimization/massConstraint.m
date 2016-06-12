%% Mass constraint
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Computes the mass constraint _cm_ and its gradient _dcmdx_ of 
% configuration _x_. A relaxation factor can be used (default value = 1).
% 
% The constraint is an inequality constraint. It is satifsfied if $c_m \geq
% 0$, that is if the mass of the structure is inferior to the maximum
% allowed mass.

function [cm,dcmdx] = massConstraint(x,relaxation)

    global mmCon
    
    if nargin < 2
        relaxation = 1;
    end
    
    cm = 0;
    dcmdx = zeros(size(x));
    
    for i = 1 : 5 : length(x)
        cm = cm + x(i+3)*x(i+4)*mmCon.rm/4;
        dcmdx(i+3) = -x(i+4)*mmCon.rm/4;
        dcmdx(i+4) = -x(i+3)*mmCon.rm/4;
    end
    
    cm = relaxation*mmCon.mMax - cm;
    
end