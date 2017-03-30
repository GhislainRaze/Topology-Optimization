

        homogenous = false;
        nMerged = 0;
        if merCon.merge && iter < oCon.iterMax
            disp('Merging mass nodes...')
            [mnodes,nMerged] = merging(mnodes);
            disp([num2str(nMerged),' mass nodes merged']);
            mmCon.n = length(mnodes);
        end
        
        nSuppressed = 0;
        if merCon.suppressZeroWidthNodes && iter < oCon.iterMax
            mnodes = suppressZeroWidthNodes(mnodes);
            nSuppressed = mmCon.n-length(mnodes);
            disp([num2str(nSuppressed),' zero-width mass nodes suppressed']);
            mmCon.n = length(mnodes);
        end
        if merCon.suppressOutsideNodes && iter < oCon.iterMax
            mnodes = suppressOutsideNodes(mnodes);
            nSuppressed = nSuppressed + mmCon.n-length(mnodes);
            disp([num2str(mmCon.n-length(mnodes)),' outside mass nodes suppressed']);
            mmCon.n = length(mnodes);
        end
        if merCon.suppressIsolatedNodes && iter < oCon.iterMax
            mnodes = suppressIsolatedNodes(mnodes,merCon.isolationRatio);
            nSuppressed = nSuppressed + mmCon.n-length(mnodes);
            disp([num2str(mmCon.n-length(mnodes)),' isolated mass nodes suppressed']);
            mmCon.n = length(mnodes);
        end
        
        if nMerged ~= 0 || nSuppressed ~= 0
            disp(['Optimization of the new structure with ',num2str(mmCon.n),' mass nodes'])
            x0 = mnodesToVector(mnodes,distrType);
            [C0,dCdx0,u0] = objectiveFunction(x0);
            relDif = 1;
            deltaX = 1;
            dt = max(2*dt,moveLimit);
        end