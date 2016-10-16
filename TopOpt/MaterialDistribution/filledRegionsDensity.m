

function rho = filledRegionsDensity(x,regions,rf)

    rho = 0;

    for h = 1 : length(regions)
        if regions(h).type == 1                              % Rectangle
            d1 = abs(x(1)-regions(h).x0(1))/regions(h).l(1)/2;
            d2 = abs(x(2)-regions(h).x0(2))/regions(h).l(2)/2;
            
            if d1 >= 1+rf || d2 >= 1+rf
                continue
            else
                if d1 <= 1-rf
                    w1 = 1;
                else
                    xi = (1+rf-d1)/(2*rf);
                	w1 = 3*xi^2-2*xi^3;
                end
                if d2 <= 1-rf
                    w2 = 1;
                else
                    xi = (1+rf-d2)/(2*rf);
                	w2 = 3*xi^2-2*xi^3;
                end
                rho = w1*w2;
                if rho == 1
                    return
                end
            end
        elseif regions(h).type == 2                          % Circle
            d = norm(x - regions(h).x0)/regions(h).r;
            if d < 1+rf
                if d > 1-rf
                    xi = (1+rf-d)/(2*rf);
                    rho = 3*xi^2-2*xi^3;
                else
                   rho = 1;
                   return
                end
            end
             
        end
    end

end