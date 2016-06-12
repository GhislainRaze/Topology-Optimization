%% Vector to Mass nodes
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Converts the mass nodes variables vector _xmnodes_ to mass nodes data

function vectorTomnodes(xmnodes,distrType)

    global mnodes mmCon

    if distrType == 3
        nd = 5;
    else
        nd = distrType+1;
    end
    
    for i = 1:length(mnodes)
        mnodes(i).x(1) = xmnodes(nd*(i-1)+1);
        mnodes(i).x(2) = xmnodes(nd*(i-1)+2);
        if distrType >= 2
            mnodes(i).theta = xmnodes(nd*(i-1)+3);
        end
        if distrType == 3
            mnodes(i).l(1) = xmnodes(nd*(i-1)+4);
            mnodes(i).l(2) = xmnodes(nd*i);
            mnodes(i).m = mmCon.rm/4*mnodes(i).l(1)*mnodes(i).l(2);
        end
    end
    
end