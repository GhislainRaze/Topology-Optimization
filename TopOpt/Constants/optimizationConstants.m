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
% This function also defines the boundary conditions.

oCon.iterMax = 1000;
oCon.relTol = 1e-6;
oCon.xTol = 1e-2;

oCon.dg = 0.01;

oCon.trueMinimum = true;
oCon.iterMinimum = 20;