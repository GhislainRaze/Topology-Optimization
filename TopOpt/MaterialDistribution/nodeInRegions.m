
function nodeInRegions = nodeInRegions(x,regions,domains)
    
    if nargin < 3
        domains = [];
    end

    nodeInRegions = false;
    
    for h = 1 : length(regions)
        if regions(h).type == 1                              % Rectangle
            nodeInRegions = min(x > regions(h).x0 - regions(h).l/2) &&...
                min(x < regions(h).x0 + regions(h).l/2);
        elseif regions(h).type == 2                          % Circle
            nodeInRegions = norm(x - regions(h).x0) < regions(h).r;
        end
        if nodeInRegions
            return
        end
    end
    
    for d = 1 : length(domains)
        if domains(d).f(x(1),x(2)) < 0
            nodeInRegions = true;
            return
        end
    end
end