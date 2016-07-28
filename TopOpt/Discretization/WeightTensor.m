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

function [w dwdx dwdy]=WeightTensor(xj,xi,dm)

wx = zeros(1,size(xj,2));
dwdrx = wx;
wy = wx;
dwdry = wx;

rx = abs(xj(1,:)-xi(1))./dm(1,:);
ry = abs(xj(2,:)-xi(2))./dm(2,:);
drdx = sign(xj(1,:)-xi(1))./dm(1,:);
drdy = sign(xj(2,:)-xi(2))./dm(2,:);

ind = 1:size(xj,2);

indx1 = ind(rx < 0.5);
indx2 = ind(and(rx >= 0.5 , rx < 1));
indx3 = ind(rx >= 1);

indy1 = ind(ry < 0.5);
indy2 = ind(and(ry >= 0.5 , ry < 1));
indy3 = ind(ry >= 1);

wx(indx1) = 2/3-4*rx(indx1).^2+4*rx(indx1).^3;
wx(indx2) = 4/3-4*rx(indx2)+4*rx(indx2).^2-4/3*rx(indx2).^3;
wx(indx3) = 0;

dwdrx(indx1) = (-8*rx(indx1)+12*rx(indx1).^2).*drdx(indx1);
dwdrx(indx2) = (-4+8*rx(indx2)-4*rx(indx2).^2).*drdx(indx2);
dwdrx(indx3) = 0;

wy(indy1) = 2/3-4*ry(indy1).^2+4*ry(indy1).^3;
wy(indy2) = 4/3-4*ry(indy2)+4*ry(indy2).^2-4/3*ry(indy2).^3;
wy(indy3) = 0;

dwdry(indy1) = (-8*ry(indy1)+12*ry(indy1).^2).*drdy(indy1);
dwdry(indy2) = (-4+8*ry(indy2)-4*ry(indy2).^2).*drdy(indy2);
dwdry(indy3) = 0;


w=wx.*wy;
dwdx=wy.*dwdrx;
dwdy=wx.*dwdry;

