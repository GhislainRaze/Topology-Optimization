function mnodes = suppressOutsideNodes(mnodes)

    global pCon

    suppressInd = [];

    for i = 1 : length(mnodes)
        if max(mnodes(i).x <= [0;-pCon.Ly/2]) || max(mnodes(i).x >= [pCon.Lx;pCon.Ly/2]) ||...
             nodeInRegions(mnodes(i).x,[pCon.holes,pCon.filledRegions],pCon.domains)
            suppressInd = [suppressInd;i];
        end
    end
        
    mnodes(suppressInd) = [];
end