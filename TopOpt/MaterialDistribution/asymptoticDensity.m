%% Asymptotic density
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Computes the asymptotic density at coordinate _xj_ of node _xi_ with
% associated _thetai_ (rotation angle), _dmi_ (smoothing lenghts) and _mi_
% (nodal mass). The other inputs are
% 
% * _rhoMin_: the minimum density (to avoid a singular stiffness matrix)
% * _rhoMax_: the maximum density before it is penalized by the asymptotic
% density
% * _distrType_: the material distribution type
% * _computeDerivatives_: set to true if the density derivatives have to be
% comuted (default value = false)
% * _rm_: the ratio between the mass node domain dimensions and its mass,
% ie $r_m = m^I/(d_x^Id_y^I)$
%
%
% The density is first based on a kernel approximation with a cubic spline
% function.
%
% $$\rho(\mathbf{x}) = \sum_{I=1}^{n} \phi^I(\mathbf{x})m^I$$
%
%
% Densities greater than one are then penalized by the asymptotic density.
% The asymptotic density is an operator on the density that is almost
% linear up to _rhoMax_ and then decreases drastically.
%
% $$ \rho^{as}(\rho(\mathbf{x})) = \frac{a\rho(\mathbf{x})}{\rho(\mathbf{x})^b + a} $$
%
% The constants $a$ and $b$ are set so that 
%
% $$\rho^{as}(\rho_{Max}) = 1 $$ and $$\frac{d \rho^{as}}{d\rho(\mathbf{x})}(\rho_{Max}) = 0$$
%
%
% Finally a minimum density is added to avoid stiffness matrix
% singularities
%
% $$ \rho^{m}(\mathbf{x}) = \rho_{Min} +
% (1-\rho_{Min})\rho^{as}(\rho(\mathbf{x})) $$

function [rhoA,drhoAdx] = asymptoticDensity(xj,xi,thetai,dmi,rm,rhoMax,...
    distrType,computeDerivatives)

    if nargin < 7
        distrType = 1;
    end
    if nargin < 8
        computeDerivatives = false;
    end

    nn = size(xi,2);
    
    if distrType == 3
        nd = 5;
    else
        nd = distrType+1;
    end
    
    drhodx = zeros(1,nd*nn);
    
    % Density and its derivatives evaluation
    ci = cos(thetai);
    si = sin(thetai);
    cxi = (xj(1)-xi(1,:)).*ci;
    cyi = (xj(2)-xi(2,:)).*ci;
    sxi = (xj(1)-xi(1,:)).*si;
    syi = (xj(2)-xi(2,:)).*si;

    r1 = cxi + syi;
    r2 = -sxi + cyi;
        
    [wi dwidx dwidy] = WeightTensor([r1;r2],zeros(2,1),dmi);
    
    rho = 4*rm*sum(wi);
    if computeDerivatives
        drhodx(1:nd:end-nd+1) = 4*rm*(-ci.*dwidx + si.*dwidy);
        drhodx(2:nd:end-nd+2) = -4*rm*(si.*dwidx + ci.*dwidy);
        if distrType >= 2
            drhodx(3:nd:end-nd+3) = 4*rm*((-sxi+cyi).*dwidx -(cxi+syi).*dwidy);
        end
        if distrType == 3
            drhodx(4:nd:end-nd+4) = -2*rm*dwidx.*(cxi+syi)./dmi(1,:);
            drhodx(5:nd:end-nd+5) = 2*rm*dwidy.*(sxi-cyi)./dmi(2,:);
        end
    end
    
    % Asymptotic density
    b = 1/(rhoMax-1)+1;
    a = rhoMax^b/(rhoMax-1);
    
    
    rhoA = a*rho./(rho.^b+a);
    
    % Density derivatives
    if computeDerivatives
        drhoAdx = (a*(1-b)*rho.^b + a^2)./((rho.^b+a).^2).*drhodx;
    end
    

end