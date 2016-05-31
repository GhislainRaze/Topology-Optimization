%% Topology optimization
% 

clear all
close all
clc



addpath('Optimization/');
addpath('Constants/');
addpath('Discretization/');
addpath('MaterialDistribution/');
addpath('EFG/');
addpath('FEM/');
addpath('IIEFG/');
addpath('Plots/');

GlobalConst;

%% Method, material distribution and optimization algorithm
% The elastic problem can be discretized thanks to three different methods
% : the Element-Free Galerkin (EFG) method, the Finite Element Method (FEM)
% or the Improved Interpolating EFG (IIEFG) method.
%
% The material distribution 
methodChoice = 2;           % 1: EFG, 2: FEM, 3: IIEFG 
massChoice = 1;             % 1: Mass nodes, 2: Undeformable structural 
                            % members, 3: Deformable structural members 
optimChoice = 1;            % 1: Overvelde's algorithm, 2: fminunc or 
                            % fmincon, 3: Genetic algorithm
                            
                            
%% Plots and movies

plotInitial = true;         % Plots the initial mass distribution
plotFinal = true;           % Plots the final mass distribution
plotMesh = true;            % Plots the discretization mesh
makeMovie = false;          % Makes a movie of the evolution of the 
                            % material distribution (slows down the
                            % algorithms)



%% Optimization process

switch optimChoice
    
    case 1
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                  Overvelde''s algorithm               ')
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        overvelde(massChoice,methodChoice);
    case 2
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                      Matlab fminunc                   ')  
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        matlabFminunc(massChoice,methodChoice);
    case 3
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                 Matlab genetic algorithm              ')
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        matlabGa(massChoice);
end


