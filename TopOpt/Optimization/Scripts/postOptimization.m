%% Post Optimization steps
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Displays information on the optimization algorithm ending reason(s).

if relDif <= oCon.relTol
    disp(['Algorithm ended after ',num2str(iter),...
        ' iterations, maximum relative tolerance on objective function reached']);
end

if deltaX <= oCon.xTol
    disp(['Algorithm ended after ',num2str(iter),...
        ' iterations, maximum tolerance on variables reached']);
end

if iter == oCon.iterMax
    disp(['Algorithm ended after ',num2str(iter),...
        ' iterations, maximum number of iteration reached']);
end

disp(['Relative change in objective function: ', num2str(relDif)])
disp(['Maximum change in variables: ', num2str(deltaX)])