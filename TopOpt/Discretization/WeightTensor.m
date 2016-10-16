%% Weight Tensor function
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Kernel function used in the MLS shape function at coordinate _xj_ for a
% node at coordinate _xi_. _dm_ is the smoothing length.

function [w dwdx dwdy]=WeightTensor(xj,xi,dm,computeDerivatives)

wx = zeros(1,size(xj,2));
wy = wx;
dwdrx = wx;    
dwdry = wx;

rx = abs(xj(1,:)-xi(1))./dm(1,:);
ry = abs(xj(2,:)-xi(2))./dm(2,:);

if computeDerivatives
    drdx = sign(xj(1,:)-xi(1))./dm(1,:);
    drdy = sign(xj(2,:)-xi(2))./dm(2,:);
end

indx1 = false(size(xj,2),1);
indx2 = indx1;
indx3 = indx1;
indy1 = indx1;
indy2 = indx1;
indy3 = indx1;

for i = 1 : size(xj,2)
    
    if rx(i) < 0.5
        indx1(i) = true;
    elseif rx(i) < 1
        indx2(i) = true;
    else
        indx3(i) = true;
    end
    
    if ry(i) < 0.5
        indy1(i) = true;
    elseif ry(i) < 1
        indy2(i) = true;
    else
        indy3(i) = true;
    end
    
end

wx(indx1) = 2/3-4*rx(indx1).^2+4*rx(indx1).^3;
wx(indx2) = 4/3-4*rx(indx2)+4*rx(indx2).^2-4/3*rx(indx2).^3;
wx(indx3) = 0;

wy(indy1) = 2/3-4*ry(indy1).^2+4*ry(indy1).^3;
wy(indy2) = 4/3-4*ry(indy2)+4*ry(indy2).^2-4/3*ry(indy2).^3;
wy(indy3) = 0;

if computeDerivatives
    dwdrx(indx1) = (-8*rx(indx1)+12*rx(indx1).^2).*drdx(indx1);
    dwdrx(indx2) = (-4+8*rx(indx2)-4*rx(indx2).^2).*drdx(indx2);
    dwdrx(indx3) = 0;

    dwdry(indy1) = (-8*ry(indy1)+12*ry(indy1).^2).*drdy(indy1);
    dwdry(indy2) = (-4+8*ry(indy2)-4*ry(indy2).^2).*drdy(indy2);
    dwdry(indy3) = 0;
end

w=wx.*wy;
dwdx=wy.*dwdrx;
dwdy=wx.*dwdry;

