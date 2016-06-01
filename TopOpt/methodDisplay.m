%% Method display
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Displays the chosen discretization method.

function methodDisplay(choice)
switch choice
    case 1
        disp('                          EFG                          ')
    case 2
        disp('                          FEM                          ')
    case 3
        disp('                         IIEFG                         ')
end
end