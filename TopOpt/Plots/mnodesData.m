%% Mass nodes data
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Converts data from a vector _x_ to mass nodes characteristics. The inputs
% are
%
% * _nd_: the number of variables per mass node
% * _x_: the mass nodes variables vector
%
% The outputs are
%
% *_xi_: the mass nodes coordinates
% *_thetai_: the mass nodes rotation angles
% *_dmi_: the mass nodes influence domain dimensions
% *_mi_: the mass nodes masses

function [xi,thetai,dmi,mi] = mnodesData(nd,x)
    
    global mmCon;
    nmn = length(x)/nd;
    xi = zeros(2,nmn);
    thetai = zeros(1,nmn);
    dmi = xi;
    mi = thetai;

    for i = 1:nd:length(x)
        nn = ceil(i/nd);
        xi(1,nn) = x(i);
        xi(2,nn) = x(i+1);
        if nd >= 3
            thetai(nn) = x(i+2);
        else
            thetai(nn) = 0;
        end
        if nd == 5
            dmi(1,nn) = x(i+3)/2;
            dmi(2,nn) = x(i+4)/2;
            mi(nn) = dmi(1,nn)*dmi(2,nn)*mmCon.rm;
        else
            dmi(1:2,nn) = mmCon.dm;
            mi(nn) = mmCon.mi;
        end
    end 
    
end