# Topology-Optimization
Topology optimization based on the Moving Node Algorithm using EFG, FEM or IIEFG methods.

## Introduction
This set of Matlab files are used for topology optimization using the Moving Node Approach (MNA) for 2D plane strain problems. In this approach, the material distribution is decoupled from the discretization.

### Material distribution
The material distribution is used to specify where the material is. It is based on mass nodes. The density at one point can be computed thanks to a kernel approximation using cubic spline shape functions.

![](http://latex.codecogs.com/gif.latex?%24%24%5Crho%28%5Cmathbf%7Bx%7D%29%20%3D%20%5Csum_%7BI%3D1%7D%5E%7Bn%7D%20%5Cphi%5EI%28%5Cmathbf%7Bx%7D%29m%5EI%24%24)

Some corrections are added in order to avoid numerical problems.

### Discretization
The governing linear elasticity equations have to be discretized to solve numerically the problem. The discretization methods can be :
* A meshless method called Element-Free Galerkin (EFG)
* A Finite Element Method (FEM)


## Code structure

### Main file
The main file topologyOptimization.m is a Matlab script that launches the optimizers. The user can change the script to set
* The discretization method (EFG or FEM)
* The optimization variables (mass nodes, undeformable structural members or deformable structural members)
* The optimization algorithm

### Problem constants
The problem constants are defined in the `Constants\` directory. This includes
* Problem geometry
* Boundary conditions
* Material distribution constants
 
### Discretization
The discretization methods use functions from the `Discretization\`, `EFG\` and `FEM\` directories.

The files `InitEFGMesh.m` and `InitFEMMesh.m` contain the discretization constants. This includes
* Number of nodes or element per dimension
* Shape functions degree
* Number of Gauss points per cell/element
* Relative smoothing length (EFG)

The other files of these directories should not be modified directly by the user.

### Optimization
The optimizers as well as some useful functions are placed in the `Optimization\` directory. Depending on the algorithm, some specific constants can be found in this directory files.

