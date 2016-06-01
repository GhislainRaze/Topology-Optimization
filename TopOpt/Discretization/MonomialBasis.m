%% Monomial Basis
%
% Created by J.T.B. Overvelde - 18 April 2012
%
% Master's thesis - The Moving Node Approach in Topology Optimization
%
% <http://www.overvelde.com>
%
% Monomial basis at coordinate _x_ for _pn_ terms.

function [p dpdx dpdy]=MonomialBasis(x,pn)
 
if pn==1
    p=1;
    dpdx=0;
    dpdy=0;
elseif pn==2
    p=[1 x(1) x(2)]';
    dpdx=[0 1 0]';
    dpdy=[0 0 1]';
elseif pn==3
    p=[1 x(1) x(2) x(1)^2 x(2)^2 x(1)*x(2)]';
    dpdx=[0 1 0 2*x(1) 0 x(2)]';
    dpdy=[0 0 1 0 2*x(2) x(1)]';
elseif pn==4
    p=[1 x(1) x(2) x(1)^2 x(2)^2 x(1)*x(2) x(1)^3 x(2)^3 x(1)^2*x(2) x(1)*x(2)^2]';
    dpdx=[0 1 0 2*x(1) 0 x(2) 3*x(1)^2 0 2*x(1)*x(2) x(2)^2]';
    dpdy=[0 0 1 0 2*x(2) x(1) 0 3*x(2)^2 x(1)^2 2*x(1)*x(2)]';
else
    disp('chose value for pn between 1 and 4')
end
        

    