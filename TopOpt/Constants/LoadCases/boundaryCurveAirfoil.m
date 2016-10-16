%% Airfoil Boundary Curve
%
% Level-set function to describe the airfoil domain

function f = boundaryCurveAirfoil(x,y,R,dx,dy)
    
    f = R^2 - (x - dx).^2 - (y - dy).^2;
        
end