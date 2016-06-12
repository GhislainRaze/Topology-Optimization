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

function [rho,X,Y,xi] = densityField(x,distrType,n)

    global pCon mmCon;

    if nargin < 3
        n = 7;
    end
    if distrType == 3
        nd = 5;
    else
        nd = distrType+1;
    end
    
    [X,Y] = meshgrid(linspace(0,pCon.Lx,n),linspace(-pCon.Ly/2,pCon.Ly/2,n));
    dx = pCon.Lx/(n-1);
    dy = pCon.Ly/(n-1);
    rho = ones(size(X));

    
    nodelabel = 1:length(x)/nd;

    [xi,thetai,dmi,mi] = mnodesData(nd,x);

    for ii = 1 : size(X,1)
        for jj = 1 : size(X,2)
            if rho(ii,jj) > 0
                xj = [ X(ii,jj) , Y(ii,jj) ]; 

                r1 = (xi(1,:)-xj(1)).*cos(thetai)+(xi(2,:)-xj(2)).*sin(thetai);
                r2 = -(xi(1,:)-xj(1)).*sin(thetai)+(xi(2,:)-xj(2)).*cos(thetai);

                nemn=nodelabel(and(abs(r1)<dmi(1,:),abs(r2)<dmi(2,:)));

                if ~isempty(nemn)
                    xinn = xi(:,nemn);
                    thetainn = thetai(nemn);
                    dminn = dmi(:,nemn);
                    minn = mi(nemn);

                    rho(ii,jj) = asymptoticDensity(xj,xinn,thetainn,dminn,minn,...
                        mmCon.rhoMax,distrType);
                    
                else
                    rho(ii,jj) = 0;
                end
            end
        end
    end
    
    for i=1:2
        
        n = 4+(n-2)*3;
        dx = dx/3;
        dy = dy/3;
        XX = zeros(n);
        YY = XX;
        rho2 = XX;
        k=0;

        for i = 1:3:size(XX)
            k = k+1;
            l = 0;
            for j = 1:3:size(YY)
                l = l+1;
                XX(i,j) = X(k,l);
                YY(i,j) = Y(k,l);
                rho2(i,j) = rho(k,l);
                if i ~= 1
                    XX(i-1,j) = X(k,l);
                    YY(i-1,j) = Y(k,l)-dy;
                    if j ~=1
                       XX(i-1,j-1) = X(k,l)-dx;
                       YY(i-1,j-1) = Y(k,l)-dy;
                    end
                    if j~=n
                       XX(i-1,j+1) = X(k,l)+dx;
                       YY(i-1,j+1) = Y(k,l)-dy;
                    end
                end
                if i~= n
                    XX(i+1,j) = X(k,l);
                    YY(i+1,j) = Y(k,l)+dy;
                    rho2(i+1,j) = rho(k,l) + (rho(k+1,l)-rho(k,l))/3;
                    rho2(i+2,j) = rho(k,l) + 2*(rho(k+1,l)-rho(k,l))/3;
                    if j ~=1
                       XX(i+1,j-1) = X(k,l)-dx;
                       YY(i+1,j-1) = Y(k,l)+dy;
                    end
                    if j~=n
                       XX(i+1,j+1) = X(k,l)+dx;
                       YY(i+1,j+1) = Y(k,l)+dy;
                       rho2(i+1,j+1) = rho(k,l) + (rho(k+1,l+1)-rho(k,l))/3;
                       rho2(i+2,j+2) = rho(k,l) + 2*(rho(k+1,l+1)-rho(k,l))/3;
                    end
                end
                if j ~= 1
                    XX(i,j-1) = X(k,l)-dx;
                    YY(i,j-1) = Y(k,l);
                end
                if j~= n
                    XX(i,j+1) = X(k,l)+dx;
                    YY(i,j+1) = Y(k,l);
                    rho2(i,j+1) = rho(k,l) + (rho(k,l+1)-rho(k,l))/3;
                    rho2(i,j+2) = rho(k,l) + 2*(rho(k,l+1)-rho(k,l))/3;
                    if i ~=n
                       rho2(i+2,j+1) = 0.5*(rho2(i+2,j+2)+rho2(i+2,j));
                       rho2(i+1,j+2) = 0.5*(rho2(i+2,j+2)+rho2(i,j+2));
                    end
                end
            end
        end

        rho = rho2;
        X = XX;
        Y = YY;
        
        for ii = 1 : size(X,1)
            for jj = 1 : size(X,2)
                if rho(ii,jj) > 0  && ~(~mod(ii-1,3) && ~mod(jj-1,3))
                    xj = [ X(ii,jj) , Y(ii,jj) ]; 

                    r1 = (xi(1,:)-xj(1)).*cos(thetai)+(xi(2,:)-xj(2)).*sin(thetai);
                    r2 = -(xi(1,:)-xj(1)).*sin(thetai)+(xi(2,:)-xj(2)).*cos(thetai);

                    nemn=nodelabel(and(abs(r1)<dmi(1,:),abs(r2)<dmi(2,:)));

                    if ~isempty(nemn)
                        xinn = xi(:,nemn);
                        thetainn = thetai(nemn);
                        dminn = dmi(:,nemn);
                        minn = mi(nemn);

                        rho(ii,jj) = asymptoticDensity(xj,xinn,thetainn,dminn,minn,...
                            mmCon.rhoMax,distrType,false);
                    else
                        rho(ii,jj) = 0;
                    end
                end
            end
        end
    
    end
    
end