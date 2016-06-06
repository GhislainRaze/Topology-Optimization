%% Slider callback function for the density plot
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Callback function for the slider of the material distribution plot.

function sliderDensity(hObject,event,x,distrType,h)
        
    n = round(get(hObject, 'Value'));
    xn = x(:,n+1);
    figureTitle = ['Material distribution at iteration ',num2str(n)];
    
    densityPlot(xn,distrType,h,figureTitle);
    
end