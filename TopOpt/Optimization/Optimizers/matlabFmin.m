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
% distribution. The inputs are
%
% * _distrType_: the material distribution type (1, 2 or 3)
% * _method_: the discretization method type (1 or 2)



function history = matlabFmin(distrType,method)
     
    global mnodes oCon mCon mmCon pCon

    history.x = [];
    history.C = [];


    if method == 1
        [Ke,f,G,q,K]=EFGUnitMatrices();
        disp('Unit matrices computed')
        if oCon.filter && ~oCon.filterIter 
            [H,Hs] = filterInitialization(cells,mCon.nG,mmCon.rmin);
            filterEnabled = true;
            disp('Filter enabled')
        else
            H = [];
            Hs = [];
            filterEnabled = false;
        end
        objectiveFunction = @(x) complianceEFG(x,distrType,Ke,f,G,q,K,...
                H,Hs);
    elseif method == 2
        [Ke,f,ubar,K]=FEMUnitMatrices();
        disp('Unit matrices computed')
        if oCon.filter && ~oCon.filterIter 
            [H,Hs] = filterInitialization(cells,mCon.nG,mmCon.rmin);
            filterEnabled = true;
            disp('Filter enabled')
        else
            H = [];
            Hs = [];
            filterEnabled = false;
        end
        objectiveFunction = @(x) complianceFEM(x,distrType,Ke,f,ubar,K,...
                H,Hs);
    end

    x0 = mnodesToVector(mnodes,distrType);
    
    % Check mesh and mass distribution
    volFrac = mmCon.vol/pCon.vol;
    disp(['Volume fraction: ',num2str(100*volFrac),'%'])
    a = min(mmCon.dx,mmCon.dy)/max(mCon.dx,mCon.dy);
    if a < 1
        disp('Warning: the mesh is probably too coarse for the mass distribution.')
        disp(['    The ratio between mass nodes influence domain and cells/element',...
            ' size should be at least 1'])
        disp(['    Current value: ',num2str(a)])
    end

    
    tic2 = tic;
    
    if distrType <3
        opt = optimset('GradObj','on','Display','iter',...
            'MaxIter',oCon.iterMax,'OutputFcn',@outfun,...
            'TolFun',oCon.relTol,'TolX',oCon.xTol);
        fminunc(objectiveFunction,x0,opt);
    else
        lb = -inf(size(x0));
        ub = inf(size(x0));
        lb(4:5:end-1) = 0;
        lb(5:5:end) = 0;
        opt = optimset('GradObj','on','GradConstr','on','Display','iter',...
            'MaxIter',oCon.iterMax,'OutputFcn',@outfun,...
            'TolFun',oCon.relTol,'TolX',oCon.xTol);
        fmincon(objectiveFunction,x0,[],[],[],[],lb,ub,...
            @matlabMassConstraint,opt);
    end
    
    time = toc(tic2);
    hh = floor(time/3600);
    mm = floor((time-3600*hh)/60);
    ss = round(time - 3600*hh - 60*mm);
    titer = time/length(history.C);

    disp(['Total time elapsed: ',num2str(hh),'h ',num2str(mm),'m ',num2str(ss),'s '])
    disp(['Average time per iteration: ',num2str(titer),'s'])
    
    function stop = outfun(x,optimValues,state)
        stop = false;
 
        switch state
            case 'iter'
                % Concatenate current point and objective function
                % value with history.
                history.C = [history.C, optimValues.fval];
                history.x = [history.x, x];
            otherwise
        end
    end

end
 