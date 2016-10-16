%% Matlab Genetic Algorithm
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Calls Matlab's _ga_ (genetic algorithm) to optimize the material
% distribution. The inputs are
%
% * _distrType_: the material distribution type (1, 2 or 3)
% * _method_: the discretization method type (1 or 2)



function history = matlabGa(distrType,method)
    global mnodes oCon mCon mmCon pCon cells

    history.x = [];
    history.C = [];
    
    
    if method == 1
        [Ke,f,G,q]=EFGUnitMatrices();
        disp('Unit matrices computed')
        if oCon.filter && ~oCon.filterIter 
            [H,Hs] = filterInitialization(cells,mCon.nG,oCon.rmin);
            filterEnabled = true;
            disp('Filter enabled')
        else
            H = [];
            Hs = [];
            filterEnabled = false;
        end
        objectiveFunction = @(x) complianceEFG(x,distrType,Ke,f,G,q,...
                H,Hs,false,false);
    elseif method == 2
        [Ke,f,ubar]=FEMUnitMatrices();
        disp('Unit matrices computed')
        if oCon.filter && ~oCon.filterIter 
            [H,Hs] = filterInitialization(cells,mCon.nG,oCon.rmin);
            filterEnabled = true;
            disp('Filter enabled')
        else
            H = [];
            Hs = [];
            filterEnabled = false;
        end
        objectiveFunction = @(x) complianceFEM(x,distrType,Ke,f,ubar,...
                H,Hs,false,false);
    end

    % Check mesh and mass distribution
    disp(['Volume fraction: ',num2str(100*mmCon.volFrac),'%'])
    a = min(mmCon.dx,mmCon.dy)/max(mCon.dx,mCon.dy);
    if a < 1
        disp('Warning: the mesh is probably too coarse for the mass distribution.')
        disp(['    The ratio between mass nodes influence domain and cells/element',...
            ' size should be at least 1'])
        disp(['    Current value: ',num2str(a)])
    end
    
    if distrType == 3
        nd = 5;
    else
        nd = distrType+1;
    end
    
    LB = zeros(nd*length(mnodes),1);
    UB = LB;
    
    UB(1:nd:end-nd+1) = pCon.Lx;
    LB(2:nd:end-nd+2) = -pCon.Ly/2;
    UB(2:nd:end-nd+2) = pCon.Ly/2;
    if distrType >= 2
        LB(3:nd:end-nd+3) = -pi/2;
        UB(3:nd:end-nd+3) = pi/2;
    end
    if distrType == 3
        UB(4:nd:end-1) = inf;
        UB(5:nd:end) = inf;
    end
    
    x0 = mnodesToVector(mnodes,distrType);
    pop = 2*length(x0);
    
    opt = gaoptimset('Generations',oCon.iterMax,...
            'Display','iter','OutputFcn',@outfun,...
            'InitialPopulation',x0','PopulationSize',pop);
    
    tic2 = tic;
    
    if distrType <3
        ga(objectiveFunction,nd*length(mnodes),[],[],[],[],LB,UB,[],opt);
    else
        ga(objectiveFunction,nd*length(mnodes),[],[],[],[],LB,UB,...
            @matlabMassConstraint,opt);
    end
    
    
    time = toc(tic2);
    hh = floor(time/3600);
    mm = floor((time-3600*hh)/60);
    ss = round(time - 3600*hh - 60*mm);
    titer = time/length(history.C);
    
    disp(['Total time elapsed: ',num2str(hh),'h ',num2str(mm),'m ',num2str(ss),'s '])
    disp(['Average time per iteration: ',num2str(titer),'s'])
    
    function [state,options,optchanged] = outfun(options,state,flag,interval)
        state.StopFlag = '';
        optchanged = false;
        switch flag
            case 'iter'
                [C,Imax] = max(state.Score);
                % Concatenate current point and objective function
                % value with history.
                history.C = [history.C, C];
                history.x = [history.x, state.Population(Imax,:)'];
            otherwise
        end
    end
end
