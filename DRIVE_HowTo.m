% Tutorial script for how to set up FE problem with the present toolbox.
%
%   Note that this script is not runable, as it is just a collection of
%   code snippets.
return;
%
% Available/predefined problems:
% (see respective DRIVE_ files)
%   Poisson-eq for a point source in 2.5D         -> App_DC.DRIVE_DC
%   Poisson-eq for a point source in 2D           -> Test.Drive_Poisson
%
% Available rhs:
% (see section "Define problem specific nodal geometry.")
%   Homogenous rhs                                -> 'none'
%   Analytic reference solution                   -> 'reference'
%   Point source (Dirac)                          -> 'point_exact'
%   Approx. point source by Gaussian distribution -> 'point_approx'
%       -> Don't use the last variant, as it is not properly tested yet!
%
% Available mesh-creation-routines:
% (see section "Set up mesh.")
%   Create a rectangle splitted in triangles
%       -> own routine: Mesh.createUnitCubeMesh.m
%   Create a rectangle splitted in triangles including a rhomub structure
%       -> own routine: Mesh.createRhombMesh.m
%   Create a custom mesh using Gmsh, including TX/RX points and topography
%       -> requires external software Gmsh
%       -> TX/RX coordinates and topography info (x,y-coords) have to be
%          set up in advance
%
% Available boundary conditions:
% (see section "Define problem specific domain information.")
%   homogenous/inhomogenous Dirichlet             -> 'dirichlet'
%   homogenous/inhomogenous Neumann               -> 'neumann'
%
% General remarks:
%   Most of the usual required steps to set up the FE forward problem are
%   summarized within the <package>.init(...) and <package>.assemble(...) 
%   routines.
%   If possible, it is recommended to use these (they are fairly tested).
%
%   To be more flexible in defining your own problem you may need to 
%   exclude specific steps.
%
% Coordinate system (Kartesian):
% Note that topography therefore typically has a negative y-component.
%  0,0 ------> 
%      |      x
%      |
%      |
%      v
%       y

%% Set up script.

% Clean up workspace and set verbosity.
clean();                            % clears and closes everything
warning('on');
verbosity = pick(2, false, true);   % each function supporting verbosity 
                                    %  will print out some statistics

%% Define problem specific nodal geometry.

% Define sources (there are several source types implemented):
% 1) homogenous rhs
TX = struct();
TX.type = 'none';                   % the type field is always required!

% 2.1) singel Dirac point source
TX.type = 'point_exact';
TX.coo = [0, 0];                    % location [x, y]
TX.val = 1;                         % source strengh
% 2.2) multiple Dirac point sources
TX.type = 'point_exact';
TX.coo = [0, 0; 1, 1];              % location [x_1, y_1; x_2, y_2; ...]
TX.val = [1; -1];

% 3) an analytic reference solution
%   -> this source type was introduced to verify the code, as it allows for
%   convergence studies.
TX.type = 'reference';
TX.ref_sol_u = RefSol.getSin(); 
    % Analytic solution for the current problem
    % -> required for setting up exact Dirichlet/Neumann boundary
    %  conditions and for comparing the FE solution to that
    % -> !ONLY! use the solutions provided in the RefSol package
    %  as they all provide a consistent output
    % If you want to add own functions, place them at the package folder
    %  and mimic the available functions!
    % Note that some functions may require additional input arguments
    %  (e.g. coordinates; RefSol.getPoisson2D(TX.coo))
TX.ref_sol.f = @(X, Y) TX.ref_sol_u.L(X, Y) + TX.ref_sol_u.f(X, Y);
    % Analytic solution of the current forward problem 
    %  (e.g. in the above example: -\nabla²(u) + u)
    % -> This function handle will be used to set up the rhs of the FE
    %  problem

% Define receiver locations:
% They may also be included as nodes in the mesh ('create_gmsh') and are 
% required to set up the FE interpolation (see FeK.getInterpolation)
RX = struct();
RX.coo = [-0.5, 0; -0.25, 0];         % location [x_1, y_1; x_2, y_2; ...]
% If no observation is required just define an empty vecror.
RX.coo = [];
    
% Define topography information:
% x-y location of known points from surface
topo = struct();
topo.coo = [TX.coo; RX.coo; -3, 1; 3, 2]; % location [x_1, y_1; x_2, y_2; ...]

%% Set up mesh.

% There are several ways to set up an appropriate mesh for your problem:
% 1) create a simple mesh by your own:
mesh_type = pick(2, 'rhomb', 'cube');  % two available rectangular domains, 
                                       %  splitted up in triangles
                                       % see +Mesh.createRhombMesh.m
                                       %     +Mesh.createUnitCubeMesh.m
ref_steps = 1;                         % numer of uniform mesh refinemets
%             xmin xmax ymin ymax 
domain_bnd = [-1,  1,   -1,  1];       % domain boundaries
mesh = Mesh.initMesh(mesh_type, 'bnd', domain_bnd, ...
        'ref', ref_steps, 'verbosity', verbosity);
    
% 2) create a custom mesh with Gmsh (used for the DC application):
% Note that given TX/RX positions and topography points are assumed to form
%  the surface (i.e. the boundary at ymin). Hence, only if any TX/RX is
%  horizontally (i.e. w.r.t. x-coord) aligned with a topography point, this
%  TX/RX can be placed at the domain interior!
% Note that the domain boundary may be extended by several kilometers to
%  ensure that homogenous Dirichlet BC are approximately correct.
mesh_type = 'gmsh_create';
mesh = Mesh.initMesh(mesh_type, 'bnd', domain_bnd, ...
    'ref', ref_steps, 'verbosity', verbosity, ...
    'topo', topo.coo, 'TX', TX.coo, ...
    'RX', RX.coo, 'dom_name', '<name-of-halfspace>');

% 3) load an external mesh (currently only Gmsh is supported)
mesh_type = 'gmsh_load';
msh_name = '<meshname>.msh';
mesh = Mesh.initMesh('gmsh_load', 'name', msh_name);

%% Define problem specific domain information.

% Note that for external meshes some of the following information may only
% be obtained by first loading the mesh (and check the info contained)
%  mesh = Mesh.loadGmsh('<filename.msh>');

% Set up boundary conditions:
bnd = struct();
bnd.type = {'dirichlet', 'neumann'};        % boundary type identifier
    % For different parts of the boundary, specific boundary conditions can
    % be set.
    % In this cell struct only the different !types! are listed.
bnd.name = {'ymin', 'ymax', 'xmin', 'xmax'};% boundary part identifier
    % This cell struct contains the boundary part names (as they result e.g.
    % from the construction by the own mesh creation routines or as they 
    % need to be contained in the external mesh).
%         ymin ymax  xmin xmax 
bnd.val = {{[];   0;  0;   0}, ...   % 1 for the Dirichlet part
           {0;  [];   [];  []}}; ... % 2 for the Neumann   part
    % After setting up the boundary type and part identifier, the 
    % respective value are asigned.
    %   The outermost cell contains as many elements as bnd.type are set
    %   The innermost cell contains as many elements as bnd.name are set
    % Note that '[]' has to be set at every boundary part which is already
    % covered by an other boundary part.
    %
    % Note that a value might also be a function handle to an analytic
    % solution 
%           xmin            xmax            ymin            ymax     
bnd.val = {{TX.ref_sol_u.J;             [];             []; TX.ref_sol_u.J}, ...
           {[];             TX.ref_sol_u.f; TX.ref_sol_u.f;             []}};
    % (e.g. see Test.DRIVE_FEvsRefSol)
bnd.quad_ord = 1;                    % set quadrature order, required for 
                                     % inhomogenous Neumann BC
  
% Define domain (i.e. areas of const. conductivity) names (identifiers) and 
% the respective parameter value:
info = struct();
info.name =  {'background', 'anomaly'};     % domain part identifier
info.val = 1./[1e3,          1e-3];         % domain part parameter value
param = Param.initParam(mesh, info);        
    % A parameter vector will be created.
    % Also consistency is checked - i.e. if number of domains and its
    %  identifier are existing.
    % Note, for the own meshes 'cube' and 'rhomb' the only parameter name
    % is defined as 'entire'.

%% Set up FEM.

% Set order of Lagrange elements.
FE_order = pick(2, 1, 2);

% Summarize parameter.
fwd_params = struct();
fwd_params.TX = TX;
fwd_params.RX = RX;
fwd_params.bnd = bnd;
fwd_params.FE_order = FE_order;
  
%% Assemble FE problem.

% 1) use the predefined problem assembling subroutines:
% 2D DC.
[fe, sol] = App_DC.assembleDC2D(mesh, param, fwd_params, verbosity);

% 2.5D DC.
% Define type of numerical integration approach.
fwd_params.FT_type = pick(3, 'Boerner', 'Bing', 'Xu');
[fe, sol, FT_info] = App_DC.assembleDC25D(mesh, param, fwd_params, verbosity);

% 2) assemble the problem by your own, which requires the following steps:
% Set up FE structure.
fe = FeL.initFiniteElement(order, mesh, RX, verbosity);
% Set up BC.
bnd = FeL.assignBC(bnd, fe, mesh, param);
% Set up system matrix (here Laplace).
sol = struct();
sol.A = FeL.assembleStiff(fe, mesh, param, verbosity);
% Set up rhs vector.
sol.b = FeL.assembleRHS(fe, mesh, TX, verbosity);
% Handle boundary conditions.
sol = FeL.treatBC(fe, mesh, sol, bnd, verbosity);

%% Solve FE problem.

% 2D DC.
u_FE = FeL.solveFwd(sol, fe, verbosity);

% 2.5D DC.
u_FE = App_DC.solveDC25D(fe, sol, FT_info, verbosity);

%% Obtain solution at RX positions.

u_RX = fe.I * u_FE;

%% Visualize solution.

Plot.plotGradient(fe, mesh, u_FE, 'sign', 'pos', 'param', param);
title('Solution gradient with parameter model.');

Plot.plotSolution(fe, mesh, u_FE, 'style', '3D');
title('Solution.');

Plot.plotMesh(mesh, param);            % mainly used for debugging