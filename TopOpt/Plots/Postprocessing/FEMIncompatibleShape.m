function phi = FEM1DShape(x,nn)

    phi = zeros(1,nn);

    if nn == 2
        phi(1) = 1/2*(1-x);
        phi(2) = 1/2*(1+x);
    elseif nn == 3
        phi(1) = 1/2*(x^2-x);
        phi(2) = (1-x^2);
        phi(3) = 1/2*(x^2+x);
    elseif nn == 4
        phi(1) = 1/16*(1-x)*(9*x^2-1);
        phi(2) = 9/16*(1-x^2)*(1-3*x);
        phi(3) = 9/16*(1-x^2)*(1+3*x);
        phi(4) = 1/16*(1+x)*(9*x^2-1);
    else
       disp('Invalid border element') 
    end
    
end
