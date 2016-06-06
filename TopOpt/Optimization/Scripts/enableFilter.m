%% Enable Filter
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Enables the filter if the criteria are met (see optimizationConstants).
% The convolution matrix is computed and the objective function is
% modified.


if oCon.filter && (iter == oCon.filterIter || abs(relDif) < oCon.relTolFilter) || ...
   oCon.filter && abs(relDif) < oCon.relTol && ~filterEnabled 
    
    [H,Hs] = filterInitialization(cells,mCon.nG,mmCon.rmin);
    if method == 1
        objfun = @(x) complianceEFG(x,distrType,Ke,f,G,q,H,Hs);
    elseif method == 2
        objfun = @(x) complianceFEM(x,distrType,Ke,f,ubar,H,Hs);
    end
    
    if abs(relDif) < oCon.relTol
        relDif = 1;
    end
    
    filterEnabled = true;
    disp('Filter enabled')
end