%% Overvelde's Optimization Algorithm
%
% Calls Overvelde's algorithm to optimize the material distribution.
%
% Overvelde's algorithm is a flow-like optimization algorithm based on a
% steepest descend and a fixed step.
%
% The flow-like equation
%
% $$ \frac{\partial^2 x_i^I(t)}{\partial t^2} + c \frac{\partial
% x_i^I(t)}{\partial t} = -\frac{\partial C(x_i^I(t))}{\partial x_i^I(t)} $$
%
% is discretized in the following way (Euler explicit time integration)
%
% $$ v_i^I(t^{k+1}) = v_i^I(t^k) - \Delta t \left(c v_i^I(t^k) + \frac{\partial 
% C(x_i^I(t^k))}{\partial x_i^I(t^k)}\right)$$
%
% $$ x_i^I(t^{k+1}) = x_i^I(t^k) + \Delta t v_i^I(t^{k+1})$$


function overvelde(distrType,method)
GlobalConst

if distrType == 3
    nd = 5;
else
    nd = distrType+1;
end

c = 0.8;
iterMax = 100;
relDif = 1;
relDifMax = 1e-6;
iter = 0;
if method == 1
    [Ke,f,G,q]=EFGUnitMatrices();
elseif method ==2
    [Ke,f,ubar]=FEMUnitMatrices();
end
disp('Unit matrices computed')

dCdx = zeros(nd*length(mnodes),1);
v = dCdx;
dt = pCon.Lx/(mCon.nx-1);
h = figure;
while relDif > relDifMax && iter <= iterMax
    if method == 1
        [~,C(iter+1),dCdx]=EFG(Ke,f,G,q,distrType);
    elseif method == 2
        [~,C(iter+1),dCdx]=FEM(Ke,f,ubar,distrType);
    end
    disp(['Iteration ',num2str(iter),' : Compliance = ',...
        num2str(C(iter+1)), ' ; Gradient norm = ',num2str(norm(dCdx))])

    if ~mod(iter,5)
        h = densityPlot(h);
    end
    
    for i = 1:nd:length(dCdx)
        nm = ceil(i/nd);
        
        normC = sqrt(dCdx(i)^2 + dCdx(i+1)^2);
        if normC == 0
            normC = 1;
        end
        dCs = dCdx(i:i+1)/normC;

        v(i) = v(i)-dt*(c*v(i)+dCs(1));
        v(i+1) = v(i+1)-dt*(c*v(i+1)+dCs(2));
        mnodes(nm).x(1) = mnodes(nm).x(1) + v(i)*dt;
        mnodes(nm).x(2) = mnodes(nm).x(2) + v(i+1)*dt;
        if distrType >= 2
            v(i+2) = (1-c)*v(i+2)-dt*dCs(3);
            mnodes(nm).theta = mnodes(nm).theta + v(i+2)*dt;
        end
        if distrType == 3
            v(i+3) = (1-c)*v(i+3)-dt*dCs(4);
            v(i+4) = (1-c)*v(i+4)-dt*dCs(5);
            mnodes(nm).l(1) = mnodes(nm).l(1) + v(i+3)*dt;
            mnodes(nm).l(2) = mnodes(nm).l(2) + v(i+4)*dt;
        end
        
        
    end
    if iter ~= 0
        relDif = abs((C(iter+1)-C(iter))/C(iter+1));
    end
    iter = iter + 1;
end
% 
densityPlot;
figure
plot(0:iter-1,C)
set(gca,'fontsize',20)
grid on
xlabel('Iteration (-)')
ylabel('Compliance')
title('Compliance convergence')
end