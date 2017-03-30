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
        objectiveFunction = @(x) complianceEFG(x,distrType,Ke,f,G,q,H,Hs);
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
        objectiveFunction = @(x) complianceFEM(x,distrType,Ke,f,ubar,H,Hs);
    end
    
    x0f = [];
    x0 = mnodesToVector(mnodes,distrType);
    
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
        Lmin = min(2*mCon.dx,2*mCon.dy);
        LB(4:nd:end-1) = Lmin;
        UB(4:nd:end-1) = inf;
        LB(5:nd:end) = Lmin;
        UB(5:nd:end) = inf;
    end
    
    tic2 = tic;
    
    if distrType <3
        % Setting 'LargeScale' to 'off' seems to be faster and to yield
        % lower compliance values (but there are more iterations)
        opt = optimset('GradObj','on','Display','iter',...
            'OutputFcn',@outfun,'LargeScale','off');
        x0 = fminunc(objectiveFunction,x0,opt);
        if filterEnabled && oCon.filterIter
            oCon.iterMax = oCon.iterMax - length(history.C);
            fminunc(objectiveFunction,x0,opt);
        end
    else
        % The algorithms 'interior-point' and 'sqp' work fine,
        % 'interior-point' is faster
        opt = optimset('GradObj','on','GradConstr','on',...
            'TolCon',oCon.relaxation-1,'Display','iter',...
            'OutputFcn',@outfun,'Algorithm','interior-point');
        x0 = fmincon(objectiveFunction,x0,[],[],[],[],LB,UB,...
            @matlabMassConstraint,opt);
        if filterEnabled && oCon.filterIter
            oCon.iterMax = oCon.iterMax - length(history.C);
            fmincon(objectiveFunction,x0,[],[],[],[],LB,UB,...
                @matlabMassConstraint,opt);
        end
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
                history.x = [history.x, [x ; x0f]];
                
                if ~filterEnabled && optimValues.iteration == oCon.filterIter && oCon.filter
                    [H,Hs] = filterInitialization(cells,mCon.nG,oCon.rmin);
                    filterEnabled = true;
                    disp('Filter enabled')
                    if method == 1
                        objectiveFunction = @(x) complianceEFG(x,distrType,Ke,f,G,q,...
                            H,Hs,false);
                    elseif method == 2
                       objectiveFunction = @(x) complianceFEM(x,distrType,Ke,f,ubar,...
                            H,Hs,false); 
                    end
                    stop = true;
                end
                
                if optimValues.iteration > 1 && abs(history.C(end-1)-history.C(end))/history.C(end-1) < oCon.relTol
                     if ~filterEnabled && oCon.filter
                        [H,Hs] = filterInitialization(cells,mCon.nG,oCon.rmin);
                        filterEnabled = true;
                        disp('Filter enabled')
                        if method == 1
                            objectiveFunction = @(x) complianceEFG(x,distrType,Ke,f,G,q,...
                                H,Hs,false);
                        elseif method == 2
                           objectiveFunction = @(x) complianceFEM(x,distrType,Ke,f,ubar,...
                                H,Hs,false); 
                        end
                        stop = true;
                     else
                        stop = true;
                        disp(['Algorithm ended after ',num2str(optimValues.iteration),...
                        ' iterations, maximum relative tolerance on objective function reached']);
                     end
                    
                end    
                if optimValues.iteration == oCon.iterMax     
                    stop = true;
                    disp(['Algorithm ended after ',num2str(optimValues.iteration),...
                    ' iterations, maximum number of iterations reached']);
                end
                if optimValues.iteration > max(1,oCon.iterMax/20) && max(abs(history.x(:,end-1)-x)) < oCon.xTol
                    stop = true;
                    disp(['Algorithm ended after ',num2str(optimValues.iteration),...
                    ' iterations, maximum tolerance on variables reached']);
                end
            otherwise
        end
    end

end
 