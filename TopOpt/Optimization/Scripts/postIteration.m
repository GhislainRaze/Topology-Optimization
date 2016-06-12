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



if oCon.continuation
    oCon.p = continuation(oCon.p,iter,oCon.pMax);
end

if distrType < 3
    disp(['Iteration ',num2str(iter),' : Compliance = ',...
    num2str(C0p), ' ; Gradient norm = ',num2str(norm(dCdx0p))])
else
    cm = massConstraint(x0,oCon.relaxation);
    C0p = C0p + oCon.mu*log(cm);
    cm = massConstraint(x0);
    cm = cm - max(0,cm);
    if cm > 0
        oCon.mu = 2*oCon.mu;
    else
        oCon.mu = 0.9*oCon.mu;
    end
    disp(['Iteration ',num2str(iter),' : Compliance = ',...
    num2str(C0p), ' ; Gradient norm = ',num2str(norm(dCdx0p)),...
    ' ; Constraint violation = ',num2str(abs(cm))])
end


dCdx0 = dCdx0p;
relDif = abs((C0-C0p)/C0p);
C0 = C0p;
deltaX = max(abs(x0-x0p));
x0 = x0p;
time = time + toc(tic2);
history.C = [history.C, C0];
history.x = [history.x, x0];
history.u = [history.u, u0];
history.m = [history.m, mTot];
%oCon.dg = 0.95*oCon.dg;

enableFilter;
tic2 = tic;