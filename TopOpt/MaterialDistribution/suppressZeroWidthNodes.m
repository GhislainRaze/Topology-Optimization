

function mnodes = suppressZeroWidthNodes(mnodes)

    suppressInd = [];

    for i = 1 : length(mnodes)
        if min(mnodes(i).l) <= 0
            suppressInd = [suppressInd;i];
        end
    end
        
    mnodes(suppressInd) = [];
end