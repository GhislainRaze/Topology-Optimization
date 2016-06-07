%% Post Iteration steps
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Post iteration script. Updates and saves the computed values. Displays
% information.


dCdx0 = dCdx0p;
relDif = abs((C0-C0p)/C0p);
C0 = C0p;
deltaX = max(abs(x0-x0p));
x0 = x0p;
history.C = [history.C, C0];
history.x = [history.x, x0];
history.u = [history.u, u0];
history.m = [history.m, mTot];

if oCon.continuation
    oCon.p = continuation(oCon.p,iter,oCon.pMax);
end
disp(['Iteration ',num2str(iter),' : Compliance = ',...
    num2str(C0), ' ; Gradient norm = ',num2str(norm(dCdx0))])
enableFilter;