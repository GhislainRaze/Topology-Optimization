%% Matlab Fmin
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Calls Matlab's _fminunc_ (for mass nodes (_distrType_ = 1) or
% undeformable structural members (_distrType_ = 2)) or _fmincon_ (for
% deformable structural members (_distrType_ = 3)) to optimize the material
% distribution.


function matlabFmin(distrType,method)
% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];
searchdir = [];
 


if distrType == 3
    nd = 5;
else
    nd = distrType+1;
end


% call optimization
GlobalConst;
if method == 1
    [Ke,f,G,q]=EFGUnitMatrices();
elseif method == 2
    [Ke,f,ubar]=FEMUnitMatrices();
end
disp('Unit matrices computed')


x0 = zeros(nd*length(mnodes),1);
for i = 1 : length(mnodes)
    x0(nd*(i-1)+1) = mnodes(i).x(1);
    x0(nd*(i-1)+2) = mnodes(i).x(2); 
    if distrType >= 2
        x0(nd*(i-1)+3) = mnodes(i).theta;
    end
    if distrType == 3
        x0(nd*(i-1)+4) = mnodes(i).l(1);
        x0(nd*i) = mnodes(i).l(2);
    end
end
disp('Material distribution initialized')

opt = optimset('GradObj','on','Display','iter',...
    'MaxIter',50,'DerivativeCheck','on');
if method == 1
    objfun = @(x) complianceEFG(x,distrType,Ke,f,G,q);
elseif method == 2
    objfun = @(x) complianceFEM(x,distrType,Ke,f,ubar);
end
if distrType <3
    [x fval exitflag output] = fminunc(objfun,x0,opt);
else
    [x fval exitflag output] = fmincon(objfun,x0,opt);
end
    
end
 