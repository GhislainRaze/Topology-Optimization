%% Density field
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Computes the density field at a given configuration. The inputs are
%
% * _x_: the mass nodes variables vector
% *_distrType_: the mass nodes type
% *_n_ the starting sampling size (optional, default value = 7)
%
% The space is discretized into _n x n_ points and the densities are
% evaluated. The discretization is then refined and the densities values
% are interpolated from the coarser mesh. If the density given by
% interpolation is lower or equal to the minimum density, the density is
% not reevaluated at that point. This recursive approach aims to gain a
% little time for this 'expansive' computation.
%
% The outputs are
%
% *_rho_: the density field
% *_X_ and _Y_: the corresponding coordinates
% *_xi_: the mass nodes coordinates

function [rho,X,Y,xi] = densityField(x,distrType)

    global pCon mCon mmCon


    if distrType == 3
        nd = 5;
    else
        nd = distrType+1;
    end
    [X,Y] = meshgrid(linspace(0,pCon.Lx,mCon.nG*mCon.mx),...
        linspace(-pCon.Ly/2,pCon.Ly/2,mCon.nG*mCon.my));
    rho = zeros(size(X));

    
    nodelabel = 1:length(x)/nd;

    [xi,thetai,dmi] = mnodesData(nd,x);

    for ii = 1 : size(X,1)
        for jj = 1 : size(X,2)
                xj = [ X(ii,jj) ; Y(ii,jj) ]; 

                r1 = (xi(1,:)-xj(1)).*cos(thetai)+(xi(2,:)-xj(2)).*sin(thetai);
                r2 = -(xi(1,:)-xj(1)).*sin(thetai)+(xi(2,:)-xj(2)).*cos(thetai);

                nemn=nodelabel(and(abs(r1)<dmi(1,:),abs(r2)<dmi(2,:)));
                
                tmp = [];
                if ~isempty(nemn)
                    tmp.x = xi(:,nemn);
                    tmp.theta = thetai(nemn);
                    tmp.l = 2*dmi(:,nemn);
                end

                rho(ii,jj) = asymptoticDensity(xj,tmp,pCon.filledRegions,mmCon.rf,mmCon.rm,...
                    mmCon.rhoMax,distrType);
        end
    end
    
end