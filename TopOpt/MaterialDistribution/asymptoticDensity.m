

function [rhoA,drhoAdx] = asymptoticDensity(xj,xi,thetai,dmi,mi,rhoMin,...
    rhoMax,distrType)

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
        
        drhodx(nd*(mm-1)+1) = mi(mm)*(1-rhoMin)*(-ci*dwidx +...
                                             si*dwidy);
        drhodx(nd*(mm-1)+2) = -mi(mm)*(1-rhoMin)*(si*dwidx +...
                                             ci*dwidy);
        if distrType >= 2
            drhodx(nd*(mm-1)+3) = mi(mm)*(1-rhoMin)*...
                ((-sxi+cyi)*dwidx -(cxi+syi)*dwidy);
        end
        if distrType == 3
            mui = mi(mm)/(4*dmi(1,mm)*dmi(2,mm));
            drhodx(nd*(mm-1)+4) = (1-rhoMin)*(2*mui*dmi(2,mm)*wi -...
                                    mi(mm)*(wi+(cxi+syi)*dwidx)/(2*dmi(1,mm)));
            drhodx(nd*mm) = (1-rhoMin)*(2*mui*dmi(1,mm)*wi -...
                                    mi(mm)*(wi+(-sxi+cyi)*dwidy)/(2*dmi(2,mm)));
        end
    end
    
    rho = rhoMin + (1-rhoMin)*rho;
    
    % Asymptotic density
    b = 1/(rhoMax-1)+1;
    a = rhoMax^b/(rhoMax-1);
    
    
    rhoA = a*rho./(rho.^b+a);
    
    
    % Asymptotic density derivatives
    drhoAdrho = (a*(1-b)*rho.^b + a^2)./((rho.^b+a).^2);
    
    drhoAdx = drhoAdrho.*drhodx;

end