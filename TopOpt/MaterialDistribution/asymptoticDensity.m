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
% (nodal mass). The other inputs are :
% 
% * _rhoMin_: the minimum density (to avoid a singular stiffness matrix)
% * _rhoMax_: the maximum density before it is penalized by the asymptotic
% density
% * _distrType_: the material distribution type
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

function [rhoA,drhoAdx] = asymptoticDensity(xj,xi,thetai,dmi,mi,rhoMin,...
    rhoMax,distrType,computeDerivatives)

    if nargin < 8
        distrType = 1;
    end
    if nargin < 9
        computeDerivatives = true;
    end

    nn = size(xi,2);
    
    if distrType == 3
        nd = 5;
    else
        nd = distrType+1;
    end
    
    drhodx = zeros(1,nd*nn);
    rho = 0;
    
    % Density and its derivatives evaluation
    for mm = 1 : nn
        ci = cos(thetai(mm));
        si = sin(thetai(mm));
        cxi = (xj(1)-xi(1,mm))*ci;
        cyi = (xj(2)-xi(2,mm))*ci;
        sxi = (xj(1)-xi(1,mm))*si;
        syi = (xj(2)-xi(2,mm))*si;
        
        r1 = cxi + syi;
        r2 = -sxi + cyi;
        
        [wi dwidx dwidy] = WeightTensor([r1;r2],zeros(2,1),dmi(:,mm));
        
        
        rho = rho+mi(mm)*wi;
        if computeDerivatives
            drhodx(nd*(mm-1)+1) = mi(mm)*(-ci*dwidx +...
                                                 si*dwidy);
            drhodx(nd*(mm-1)+2) = -mi(mm)*(si*dwidx +...
                                                 ci*dwidy);
            if distrType >= 2
                drhodx(nd*(mm-1)+3) = mi(mm)*...
                    ((-sxi+cyi)*dwidx -(cxi+syi)*dwidy);
            end
            if distrType == 3
                mui = mi(mm)/(4*dmi(1,mm)*dmi(2,mm));
                drhodx(nd*(mm-1)+4) = (2*mui*dmi(2,mm)*wi -...
                                        mi(mm)*(wi+(cxi+syi)*dwidx)/(2*dmi(1,mm)));
                drhodx(nd*mm) = (2*mui*dmi(1,mm)*wi -...
                                        mi(mm)*(wi+(-sxi+cyi)*dwidy)/(2*dmi(2,mm)));
            end
        end
    end
    
    % Asymptotic density
    b = 1/(rhoMax-1)+1;
    a = rhoMax^b/(rhoMax-1);
    
    
    rhoA = a*rho./(rho.^b+a);
    
    % Minimum density
    rhoA = rhoMin + (1-rhoMin)*rhoA;
    
    % Density derivatives
    if computeDerivatives
        drhoAdrho = (a*(1-b)*rho.^b + a^2)./((rho.^b+a).^2);
        drhoAdx = drhoAdrho.*drhodx;
        drhoAdx = (1-rhoMin)*drhoAdx;
    end
    

end