%% Finite Element Method (FEM) shape functions
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Provides the shape functions of quadrangle elements at reduced coordinate
% _x_. To identify the element, the number of nodes _nn_ is also required.
%
% If the number of nodes is inconsistent, a warning message appears and the
% shape functions are identically equal to zero.

function [phi,dphidx,dphidy]=FEMShape(x,nn)

    phi=zeros(1,nn); dphidx=phi; dphidy=phi;
    
    if nn == 4                      % Linear rectangular element
        % Corner nodes
        xi = [-1 ; -1 ; 1 ; 1];
        eta = [-1 ; 1 ; -1 ; 1];
        phi = 1/4*(1+x(1)*xi).*(1+x(2)*eta);
        dphidx = 1/4*xi.*(1+x(2)*eta);
        dphidy = 1/4*eta.*(1+x(1)*xi);
    elseif nn == 8                  % Quadratic rectangular element
        % Corner nodes
        xi = [-1 ; -1 ; 1 ; 1];
        eta = [-1 ; 1 ; -1 ; 1];
        phi([1 3 6 8]) = 1/4*(1+x(1)*xi).*(1+x(2)*eta).*...
            (x(1)*xi+x(2)*eta-1);
        dphidx([1 3 6 8]) = 1/4*xi.* ((1+x(2)*eta).*...
            (x(1)*xi+x(2)*eta-1) + (1+x(1)*xi).*(1+x(2)*eta));
        dphidy([1 3 6 8]) = 1/4*eta.* ((1+x(1)*xi).*...
            (x(1)*xi+x(2)*eta-1) + (1+x(1)*xi).*(1+x(2)*eta));
        % Mid-side nodes
        eta = [-1 ; 1];
        phi([4 5]) = 1/2*(1-x(1)^2)*(1+eta*x(2));
        dphidx([4 5]) = -x(1)*(1+eta*x(2));
        dphidy([4 5]) = 1/2*(1-x(1)^2)*eta;
        xi = [-1 ; 1];
        phi([2 7]) = 1/2*(1-x(2)^2)*(1+xi*x(1));
        dphidx([2 7]) = 1/2*(1-x(2)^2)*xi;
        dphidy([2 7]) = -x(2)*(1+xi*x(1));
    elseif nn == 12                 % Cubic rectangular element
        % Corner nodes
        xi = [-1 ; -1 ; 1 ; 1];
        eta = [-1 ; 1 ; -1 ; 1];
        phi([1 4 9 12]) = 1/32*(1+x(1)*xi).*(1+x(2)*eta).*...
            (-10+9*(x(1)^2+x(2)^2));
        dphidx([1 4 9 12]) = 1/32*(xi.*(1+x(2)*eta).*...
            (-10+9*(x(1)^2+x(2)^2))+ 18*x(1)*(1+x(1)*xi).*(1+x(2)*eta));
        dphidy([1 4 9 12]) = 1/32*(eta.*(1+x(1)*xi).*...
            (-10+9*(x(1)^2+x(2)^2))+ 18*x(2)*(1+x(1)*xi).*(1+x(2)*eta));
        % Mid-side nodes
        xi = [-1 ; -1 ; 1 ; 1];
        eta = 1/3*[-1 ; 1 ; -1 ; 1];
        phi([2 3 10 11]) = 9/32*(1+xi*x(1))*(1-x(2)^2).*(1+9*eta*x(2));
        dphidx([2 3 10 11]) = 9/32*xi*(1-x(2)^2).*(1+9*eta*x(2));
        dphidy([2 3 10 11]) = 9/32*(-(1+xi*x(1))*2*x(2).*(1+9*eta*x(2))+...
                                (1+xi*x(1))*(1-x(2)^2)*9.*eta);
        xi = 1/3*[-1 ; -1 ; 1 ; 1];
        eta = [-1 ; 1 ; -1 ; 1];
        phi([5 6 7 8]) = 9/32*(1+eta*x(2))*(1-x(1)^2).*(1+9*xi*x(1));
        dphidx([5 6 7 8]) = 9/32*(-(1+eta*x(2))*2*x(1).*(1+9*xi*x(1))+...
                                (1+eta*x(2))*(1-x(1)^2)*9.*xi);
        dphidy([5 6 7 8]) = 9/32*eta*(1-x(1)^2).*(1+9*xi*x(1));
    else
        warning('Invalid element')
    end
end