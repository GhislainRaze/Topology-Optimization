
function x0 = mnodesToVector(mnodes,distrType)


if distrType == 3
    nd = 5;
else
    nd = distrType+1;
end

x0 = zeros(nd*length(mnodes),1);
for i = 1 : length(mnodes)
    x0(nd*(i-1)+1) = mnodes(i).x(1);
    x0(nd*(i-1)+2) = mnodes(i).x(2); 
    if distrType >= 2
        x0(nd*(i-1)+3) = mnodes(i).theta;
    end
    if distrType == 3
        x0(nd*(i-1)+4) = mnodes(i).l(1);
        x0(nd*i) = mnodes(i).l(2);
    end
end

end