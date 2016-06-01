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

alphax=2/dm(1);
alphay=2/dm(2);
rx = abs(xj(1)-xi(1))/dm(1);
ry = abs(xj(2)-xi(2))/dm(2);
drdx = sign(xj(1)-xi(1))/dm(1);
drdy = sign(xj(2)-xi(2))/dm(2);
if rx>1
    wx = 0;
    dwdrx = 0;
elseif rx>0.5
    wx = alphax*(4/3-4*rx+4*rx^2-4/3*rx^3);
    dwdrx = alphax*(-4+8*rx-4*rx^2)*drdx;
else
    wx = alphax*(2/3-4*rx^2+4*rx^3);
    dwdrx = alphax*(-8*rx+12*rx^2)*drdx;
end
if ry>1
    wy = 0;
    dwdry = 0;
elseif ry>0.5
    wy = alphay*(4/3-4*ry+4*ry^2-4/3*ry^3);
    dwdry = alphay*(-4+8*ry-4*ry^2)*drdy;
else
    wy = alphay*(2/3-4*ry^2+4*ry^3);
    dwdry = alphay*(-8*ry+12*ry^2)*drdy;
end

w=wx*wy;
dwdx=wy*dwdrx;
dwdy=wx*dwdry;

