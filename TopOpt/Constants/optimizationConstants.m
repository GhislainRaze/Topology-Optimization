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
% * _oCon.relTol_: 
% * _oCon.xTol_: 
% * _oCon.dg_: 
% * _oCon.trueMinimum_: 
% * _oCon.iterMinimum_:
%

% Algorithm parameters
oCon.iterMax = 1000;                            % Maximum number of iteration
oCon.relTol = 1e-6;                             % Maximum tolerance on the relative
                                                % compliance change
oCon.xTol = 1e-4;                               % Maximum tolerance on the variables
                                                % change
oCon.dg = 0.01;                                 % Initial step size

% Linesearch parameters
oCon.trueMinimum = false;                       % Search for true minimum
oCon.iterMinimum = 50;                          % Number of iterations for the search
                                                % of a true minimum
oCon.tolMinimum = 1e-3;                         % Relative tolerance for the minimum derivative
oCon.iterWolfe = 20;                            % Number of iterations for the search
                                                % of a minimum with Wolfe criteria

% Penalization
oCon.p = 3;                                     % Intermediate density penalization
oCon.continuation = false;                      % Continuation
if oCon.continuation
    oCon.pMax = mmCon.p;
    oCon.p = 1;
end

% Filtering
oCon.filter = false;                            % Density filter
oCon.rmin = 0.3;                                % Radius under which the variations
                                                % are filtered
oCon.filterIter = 200;                          % Number of iterations after which the
                                                % filter is enabled
oCon.relTolFilter = 1e-5;                       % Maximum tolerance on the relative
                                                % compliance change after which the
                                                % filter is enabled