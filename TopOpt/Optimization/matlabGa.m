%% Matlab Fmin
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Calls Matlab's _ga_ (genetic algorithm) to optimize the material
% distribution.


function matlabGa(distrType,method)
% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];
searchdir = [];
 
% call optimization
GlobalConst
InitMesh();
[Ke,f,G,q]=EFGUnitMatrices();

x0 = zeros(2*length(mnodes),1);
for i = 1 : length(mnodes)
    x0(2*i-1) = mnodes(i).x(1);
    x0(2*i) = mnodes(i).x(2); 
    LB(2*i-1) = 0;
    LB(2*i) = -pCon.Ly/2;
    UB(2*i-1) = pCon.Lx;
    UB(2*i) = pCon.Ly/2;
end
pop = 10*length(mnodes);
opt = gaoptimset('PopulationSize',pop,'Generations',100,'UseParallel','always','Display','iter');
objfun = @(x) compliance(x,Ke,f,G,q);
[x fval exitflag output] = ga(objfun,2*length(mnodes),[],[],[],[],LB,UB,[],opt);
end
