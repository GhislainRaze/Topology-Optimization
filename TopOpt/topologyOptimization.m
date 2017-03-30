%% Topology optimization
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>

clear all
close all
clc
%% Method, material distribution and optimization algorithm
% The elastic problem can be discretized thanks to two different methods
% : the Element-Free Galerkin (EFG) method or the Finite Element Method (FEM).
%
% The material distribution 
methodChoice = 2;           % 1: EFG, 2: FEM
massChoice = 3;             % 1: Mass nodes, 2: Undeformable structural 
                            % members, 3: Deformable structural members 
optimChoice = 7;            % 1: Overvelde's algorithm, 2: steepest descent,
                            % 3: conjugated gradients, 4: quasi Newton
                            % BFGS, 5: Matlab fminunc or fmincon 
                            % (recommanded), 6: Matlab ga.

%% Load case
loadCase = 'AirfoilStiff';     % The corresponding file must exist


%% Plots

plotInitial = true;         % Plot the initial mass distribution
plotFinal = true;           % Plot the final mass distribution
plotMesh = true;            % Plot the discretization mesh
plotDEvolution = true;      % Plot the density at given iteration
plotCEvolution = true;      % Plot the elements contour at given iteration
plotCompliance = true;      % Plot the compliance evolution
plotMass = true;            % Plot the mass evolution
plotDeformed = true;        % Plot the final deformed configuration                            


%% Build path
addpath('Optimization/');
addpath('Optimization/LineSearch/');
addpath('Optimization/Optimizers/');
addpath('Optimization/Scripts/');
addpath('Constants/');
addpath('Constants/LoadCases/');
addpath('Discretization/');
addpath('MaterialDistribution/');
addpath('EFG/');
addpath('FEM/');
addpath('Plots/');
addpath('Plots/Callbacks/');
addpath('Plots/Postprocessing/');
addpath('Display/');



%% Constants initialization
GlobalConst;
problemConstants(loadCase);
optimizationConstants;
mergingConstants;


%% Optimization process

switch optimChoice
    
    case 1
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                  Overvelde''s algorithm               ')
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        history = overvelde(massChoice,methodChoice);
    case 2
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                    Steepest descent                   ')  
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        history = steepestDescent(massChoice,methodChoice);
    case 3
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                  Conjugated gradients                 ')
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        history = conjugatedGradients(massChoice,methodChoice);
    case 4
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                    Quasi-Newton BFGS                  ')
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        history = quasiNewtonBFGS(massChoice,methodChoice);
    case 5
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                  Matlab fminunc/fmincon               ')  
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        history = matlabFmin(massChoice,methodChoice);
    case 6
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                 Matlab genetic algorithm              ')
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        history = matlabGa(massChoice,methodChoice);
    case 7
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('           Gradient algorithm with maximum step         ')
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        history = gradientMax(massChoice,methodChoice);    
    case 8
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('             BFGS algorithm with maximum step          ')
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        history = BFGSMax(massChoice,methodChoice); 
        
end

% p = profile('info');
save(['TopOptResults_',date])

%% Postprocessing
endPlots;