%% Gauss Integration Points
%
% Created by J.T.B. Overvelde - 18 April 2012
%
% Master's thesis - The Moving Node Approach in Topology Optimization
%
% <http://www.overvelde.com>
%
% The position and weights of _n_ integration points in one dimension for
% normalized length.

function [t,w]=ConGauss(n)

if n==1 
    t=0;
    w=2;
elseif n==2
    t=[-sqrt(1/3),sqrt(1/3)];
    w=[1,1];
elseif n==3
    t=[-0.77459667, 0, 0.77459667];
    w=[0.55555555 0.88888889 0.55555555];
elseif n==4
    t=[-0.86113631 -0.33998104 0.33998104 0.86113631];
    w=[0.34785485 0.65214515 0.65214515 0.34785485];
elseif n==5
    t=[-0.90617985 -0.53846931 0 0.53846931 0.90617985];
    w=[0.23692689 0.47862867 0.56888889 0.47862867 0.23692689];
elseif n==6
    t=[-0.93246951 -0.66120939 -0.23861918 0.23861918 0.66120939 0.93246951];
    w=[0.17132449 0.36076157 0.46791393 0.46791393 0.36076157 0.17132449];
elseif n==7
    t=[-0.94910791 -0.74153119 -0.40584515 0.0 0.40584515 0.74153119 0.94910791];
    w=[0.12948497 0.27970539 0.38183005 0.41795918 0.38183005 0.27970539 0.12948497 ];
elseif n==8
    t=[-0.96028986 -0.79666648 -0.52553241 -0.18343464 0.18343464 0.52553241 0.79666648 0.96028986];
    w=[0.10122854 0.22238103 0.31370665 0.36268378 0.36268378 0.31370665 0.22238103 0.10122854];
else
    t=NaN;
    w=NaN;
    disp('Choose different value for number of Gauss integration points')
end