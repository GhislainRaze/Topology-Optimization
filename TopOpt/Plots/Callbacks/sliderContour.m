%% Slider callback function for the contour plot
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Callback function for the slider of the elements contour evolution plot.

function sliderContour(hObject,event,x,distrType,h)
        
    n = round(get(hObject, 'Value'));
    xn = x(:,n+1);
    figureTitle = ['Elements contours at iteration ',num2str(n)];
    
    elementsContour(xn,distrType,h,figureTitle);
    
end