


function [ll,lt] = directionalComponents(mnode,theta)
    
    % Angle between the diagonal and the median
    alpha = atan(mnode.l(2)/mnode.l(1));
    
    % Angle between the median and direction theta
    beta = theta - mnode.theta;
    
    if tan(beta) < 0
        swap = true;
    else
        swap = false;
    end
    
    beta = mod(beta,pi/2);
    
    % Computing the longitudinal length
    if beta <= alpha
        ll =  mnode.l(1)/(2*cos(beta));
    else
        ll =  mnode.l(2)/(2*sin(beta));
    end
    
    % Computing the transversal length
    if beta <= pi/2-alpha
        lt =  mnode.l(2)/(2*cos(beta));
    else
        lt =  mnode.l(1)/(2*sin(beta));
    end
    
    
    if swap
        tmp = lt;
        lt = ll;
        ll = tmp;
    end
    
end