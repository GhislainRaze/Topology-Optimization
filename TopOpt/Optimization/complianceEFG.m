%% EFG Compliance
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Function used to provide the compliance and its sensitivities based on a
% EFG method to the _fminunc_, _fmincon_ or _ga_ Matlab optimizers.
%
% The inputs are
%
% * _xmnodes_: a vector containing the optimization variables
% * _distrType_: the material distribution type
%
% Furthermore one computes once and for all the following inputs (which can
% be obtained thanks to EFGUnitMatrices)
%
% * _Ke_: the unit stiffness matrices
% * _f_: the nodal force vector
% * _G_: the lagrangian multipliers matrix
% * _q_: the second member associated with lagrangian multipliers
% * _H_: the filter convolution matrix (optional)
% * _Hs_: the sum of the filter convolution matrix lines (optional)

function [C,dCdx,u] = complianceEFG(xmnodes,distrType,Ke,f,G,q,K,H,Hs,...
    constraint,computeDerivatives)
    
    global oCon
    
    if nargin < 10
        constraint = false;
    end
    if nargin < 11
        computeDerivatives = true;
    end
    
    vectorTomnodes(xmnodes,distrType);
    
    [u,C,dCdx]=EFG(Ke,f,G,q,K,distrType,H,Hs,computeDerivatives);
    
    if constraint
        [cm,dcmdx] = massConstraint(xmnodes,oCon.relaxation);
        C = C - oCon.mu*log(cm);
        if computeDerivatives
            dCdx = dCdx - oCon.mu/cm*dcmdx;
        end
    end
    
    

end