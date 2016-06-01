%% Mass display
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Displays the chosen material distribution type

function massDisplay(choice)
switch choice
    case 1
        disp('                       Mass Nodes                      ')
    case 2
        disp('               Undeformable Structural Members         ')
    case 3
        disp('                Deformable Structural Members          ')
end
end