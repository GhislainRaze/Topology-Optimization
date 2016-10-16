%% Moving Least Squares (MLS) Shape functions
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% MLS Shape function at coordinate _y_ for nodes _x_. _dm_ is the smoothing
% length and _pn_ is the polynomial order.

function [phi,dphidx,dphidy]=MLSShape(x,y,dm,pn)
 
nn=size(x,1);
pnn=0;
for i=1:pn
    pnn=pnn+i;
end
        
A=zeros(pnn); dAdx=A; dAdy=A; 
w=zeros(1,nn); dwdx=w; dwdy=w;
phi=zeros(1,nn); dphidx=phi; dphidy=phi;
for mm=1:nn
    [wi dwidx dwidy]=WeightTensor(y,x(mm,:)',dm,true);
    p=MonomialBasis(x(mm,:),pn);
    pTp=p*p';
    A=A+wi*pTp;
    dAdx=dAdx+dwidx*pTp;
    dAdy=dAdy+dwidy*pTp;
    w(mm)=wi;
    dwdx(mm)=dwidx;
    dwdy(mm)=dwidy;
end
[p dpdx dpdy]=MonomialBasis(y,pn);

A1 = A\eye(size(A));
c=A1*p;
dcdx=A1*(-dAdx*c+dpdx);
dcdy=A1*(-dAdy*c+dpdy);

for mm=1:nn
    piT=MonomialBasis(x(mm,:),pn);
    phi(mm)=c'*piT*w(mm);
    dphidx(mm)=dcdx'*piT*w(mm)+c'*piT*dwdx(mm);
    dphidy(mm)=dcdy'*piT*w(mm)+c'*piT*dwdy(mm);
end

end
