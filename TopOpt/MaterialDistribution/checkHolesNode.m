
function nodeInHoles = checkHolesNode(x,holes)
    nodeInHoles = false;
    for h = 1 : length(holes)
        if holes(h).type == 1                              % Rectangle
            nodeInHoles = min(x > holes(h).x0 - holes(h).l/2) &&...
                min(x < holes(h).x0 + holes(h).l/2) || nodeInHoles;
        elseif holes(h).type == 2                          % Circle
            nodeInHoles = norm(x - holes(h).x0) < holes(h).r ...
                || nodeInHoles;
        end
    end
end