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
% * _H_: the filter convolution matrix (optional)
% * _Hs_: the sum of the filter convolution matrix lines (optional)

function [C,dCdx,u] = complianceFEM(xmnodes,distrType,Ke,f,ubar,...
    H,Hs,computeDerivatives)
    
    if nargin < 8
        computeDerivatives = true;
    end
    
    vectorTomnodes(xmnodes,distrType);

    [u,C,dCdx]=FEM(Ke,f,ubar,distrType,H,Hs,computeDerivatives);
    
end