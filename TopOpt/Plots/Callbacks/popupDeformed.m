%% Popup menu callback function for the deformed plot
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Callback function for the popup menu of the deformed configuration
% evolution plot.

function popupDeformed(popup,event,data,editer)
        
    n = get(popup, 'Value');
    [mag,success] = str2num(get(editer, 'String'));
    if success
        deformedPlot(data,n,mag)
    else
        warndlg('Please insert a number in the magnifying factor box.')
    end
end