%% Cubic approximation
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Interpolates a function by a cubic polynomial and returns its minimum.
% Information from two points is required.
% The inputs are
%
% * _x1_ and _x2_, the points abscissa
% * _f1 and _f2_, the function values at these points
% * _g1_ and _g2_, the function derivatives at these points
%
% If the points are the same, the function returns that point as a result.
% if the function is linear, or quadratic with a negative curvature, the 
% function returns the point _x1_ if _f1_ is smaller than _f2_, and _x2_ 
% otherwise. Finally, if the cubic function has no minimum, the function
% extends the search interval by a factor two.

function xs = cubicApproximation(x1,x2,f1,f2,g1,g2)

    if x1 == x2
        xs = x1;
    else
        
        
        b = g1;
        c = 3*(f2-f1)/(x2-x1)^2 - (2*g1+g2)/(x2-x1);
        d = 2*(f1-f2)/(x2-x1)^3 + (g1+g2)/(x2-x1)^2;
        
        if d == 0               % Quadratic function
            
            if c <= 0           % Linear function
                if f1 < f2
                    xs = x1;
                else
                    xs = x2;
                end
            else
                xs = (- b + 2*c*x1 - 3*d*x1^2)/(2*c - 6*d*x1);
            end
        elseif (c/(3*d))^2-b/(3*d) >= 0                    % Cubic function
            xs1 = x1 - c/(3*d) + sqrt((c/(3*d))^2-b/(3*d));
            xs2 = x1 - c/(3*d) - sqrt((c/(3*d))^2-b/(3*d));
            
            if 2*c + 6*d*(xs1-x1) > 0
                xs = xs1;
            else
                xs = xs2;
            end
        else
            if f2 < f1
                xs = 2*(x2-x1);
            else
                xs = -2*(x2-x1);
            end
        end
    end    
end