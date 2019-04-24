function mnodes = suppressIsolatedNodes(mnodes,isolationRadius)

    ispaired = false(length(mnodes),1);
    for i = 1 : length(mnodes)
        if ~ispaired(i)
            for j = 1 : length(mnodes)
                if j ~= i
                    xj = abs([cos(mnodes(i).theta) sin(mnodes(i).theta) ;
                        -sin(mnodes(i).theta) cos(mnodes(i).theta)]*(mnodes(j).x - mnodes(i).x));

                    if xj(1) < isolationRadius*0.5*mnodes(i).l(1) && xj(2) < isolationRadius*0.5*mnodes(i).l(2)
                        ispaired(j) = true;
                        ispaired(i) = true;
                        break
                    end
                end
            end
        end
    end
        
    mnodes(~ispaired) = [];
end