%% Problem Constants
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% This script creates the set of optimization constants _oCon_. The
% optimization constants are
%
% * _oCon.iterMax_: the maximum number of iteration
% * _oCon.relTol_: the maximum tolerance on the relative compliance change
% * _oCon.xTol_: the maximum tolerance on the variables change
% * _oCon.dg_: the initial step length
%
% * _oCon.trueMinimum_: true if the true minimum is searched
% * _oCon.iterMinimum_: the number of iterations to search the minimum
% * _oCon.tolMimimum_: the relative tolerance on the minimum derivative
% * _oCon.iterWolfe_: the iterations to find the minimum according to Wolfe
% criteria
%
% * _oCon.p_: the intermediate density penalization factor
% * _oCon.continuation_: true if the continuation strategy is used
%
% * _oCon.filter_: true if a filter is used
% * _oCon.rmin_: radius under which the mass variations are filtered
% * _oCon.filterIter_: number of iterations after which the filter is enabled
% * _oCon.relTolFilter_: maximum tolerance on the relative compliance change
% after which the filter is enabled
%
% * _oCon.relaxation_: relaxation factor for the mass constraint (only
% relevant if deformable strucutral members are used)
% * _oCon.mu_: factor for the log barrier (only relevant if deformable
% structural members are used)

% Algorithm parameters
oCon.iterMax = 500;                             % Maximum number of iteration
oCon.relTol = 1e-6;                             % Maximum tolerance on the relative
                                                % compliance change
oCon.xTol = 1e-6;                               % Maximum tolerance on the variables
                                                % change
oCon.cTol = 1e-5;
oCon.dg = min(pCon.Lx,pCon.Ly)/10;              % Initial step size

% Linesearch parameters
oCon.trueMinimum = false;                       % Search for true minimum
oCon.iterMinimum = 50;                          % Number of iterations for the search
                                                % of a true minimum
oCon.tolMinimum = 1e-3;                         % Relative tolerance for the minimum derivative
oCon.iterWolfe = 20;                            % Number of iterations for the search
                                                % of a minimum with Wolfe criteria

% Penalization
oCon.p = 1;                                     % Intermediate density penalization
oCon.continuation = false;                      % Continuation
if oCon.continuation
    oCon.pMax = oCon.p;
    oCon.p = 1;
end

% Filtering
oCon.filter = false;                            % Density filter
oCon.rmin = 0.3;                                % Radius under which the variations
                                                % are filtered
oCon.filterIter = 1000;                          % Number of iterations after which the
                                                % filter is enabled
oCon.relTolFilter = 1e-6;                       % Maximum tolerance on the relative
                                                % compliance change after which the
                                                % filter is enabled
                                                
% Mass constraint                                                
oCon.relaxation = 1.01;                        % Relaxation factor for the mass constraint
oCon.mu = 0.01;                               % Factor for the log barrier (should be small)