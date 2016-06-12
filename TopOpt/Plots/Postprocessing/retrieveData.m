%% Retrieve Data
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Retrieves displacement and total mass after a Matlab's algorithm is used.

function history = retrieveData(history,method,distrType)

    history.u = [];
    history.m = [];

    % Unit matrices
    if method == 1
        [Ke,f,G,q,K] = EFGUnitMatrices();
    elseif method == 2
        [Ke,f,ubar,K] = FEMUnitMatrices();
    end
    
    % Compute nodal displacements and total mass
    for i = 1 : length(history.C)
        
        vectorTomnodes(history.x(:,i),distrType);
        
        if method == 1
            [u,~,~,mTot] = EFG(Ke,f,G,q,K,distrType,[],[],false);
        elseif method == 2
            [u,~,~,mTot] = FEM(Ke,f,ubar,K,distrType,[],[],false);
        end
        
        history.u = [history.u,u];
        history.m = [history.m,mTot];
    end
end