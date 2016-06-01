%% FEM Compliance
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Function used to provide the compliance and its sensitivities based on a
% FEM analysis to the _fminunc_, _fmincon_ or _ga_ Matlab optimizers.
%
% The inputs are
%
% * _xmnodes_: a vector containing the optimization variables
% * _distrType_: the material distribution type
%
% Furthermore one computes once and for all the following inputs (which can
% be obtained thanks to FEMUnitMatrices)
%
% * _Ke_: the unit stiffness matrices
% * _f_: the nodal force vector
% * _ubar_: the imposed nodal displacements

function [C,dCdx] = complianceFEM(xmnodes,distrType,Ke,f,ubar)

    GlobalConst
    if distrType == 3
        nd = 5;
    else
        nd = distrType+1;
    end
    
    for i = 1:length(mnodes)
        mnodes(i).x(1) = xmnodes(nd*(i-1)+1);
        mnodes(i).x(2) = xmnodes(nd*(i-1)+2);
        if distrType >= 2
            mnodes(i).theta = xmnodes(nd*(i-1)+3);
        end
        if distrType == 3
            rm = mnodes(i).m/(mnodes(i).l(1)*mnodes(i).l(2));
            mnodes(i).l(1) = xmnodes(nd*(i-1)+4);
            mnodes(i).l(2) = xmnodes(nd*i);
            mnodes(i).m = rm*mnodes(i).l(1)*mnodes(i).l(2);
        end
    end
    [~,C,dCdx]=FEM(Ke,f,ubar,distrType);

end