%% Boundary Conditions Parameters
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% This functions takes a set of problem constants _pCon_ to compute
% the parameters associated with the boundary conditions defined
% in the <./problemConstants.html Problem Constants> file. 
%
% This includes:
%
% * Number of boundary conditions of each type
% * Length of line boundaries
% * Parametrization of each line boundary as a function of the length from
% the origin

function pCon = boundaryConditionsParameters(pCon)
    
    % Line force
    pCon.nlLoad = length(pCon.lLoad);
    for i = 1 : pCon.nlLoad
        % Line length
        pCon.lLoad(i).length = norm(pCon.lLoad(i).x(2,:)-pCon.lLoad(i).x(1,:));
        % Line parametrization
        pCon.lLoad(i).param = @(x) [interp1([0 pCon.lLoad(i).length],[pCon.lLoad(i).x(1,1) pCon.lLoad(i).x(2,1)],x,'linear');
                          interp1([0 pCon.lLoad(i).length],[pCon.lLoad(i).x(1,2) pCon.lLoad(i).x(2,2)],x,'linear')];

    end
    
    % Point force
    pCon.npLoad = length(pCon.pLoad);

    % Line essential boundary condition
    pCon.nlbc = length(pCon.lbc);
    for i = 1 : pCon.nlbc
        % Line length
        pCon.lbc(i).length = norm(pCon.lbc(i).x(2,:)-pCon.lbc(i).x(1,:));
        % 
        pCon.lbc(i).dir = [pCon.lbc(i).x(2,:)-pCon.lbc(i).x(1,:)]'/pCon.lbc(i).length;
        % Line parametrization
        pCon.lbc(i).param = @(x) [interp1([0 pCon.lbc(i).length],[pCon.lbc(i).x(1,1) pCon.lbc(i).x(2,1)],x,'linear');
                              interp1([0 pCon.lbc(i).length],[pCon.lbc(i).x(1,2) pCon.lbc(i).x(2,2)],x,'linear')];
    end
    
    % Point essential boundary condition
    pCon.npbc = length(pCon.pbc);
end