classdef testLagrange < matlab.unittest.TestCase
    % Tests Lagrange element implementation with MATLAB unittest.
    %
    % To run tests from project folder use:
    %   runtests('Test.testLagrange')

    properties (TestParameter)
        
        % Parameters for source, domain, boundary condition:
        %    test number                , ...
        %    source type     char       , ...
        %    ref. sol. field fun handle , ...
        %    x               [min, max] , ...
        %    y               [min, max] , ...
        %    bnd type        char       , ...
        %    known L2-error             , ... % known from successful runs 
        %    known H1-error             , ... % known from successful runs
        %    verbosity       locigal
        
        factory = {
            {1, 'reference',   @RefSol.getSin,       [-4, 4], [-4, 4], 'pure', ...
             [-2.1309, -2.1010; ...
              -2.5899, -2.9630], ...
             [-1.0233, -0.9733; ...
              -2.0763, -2.3718], true}, ...
            {2, 'reference',   @RefSol.getSin,       [-4, 4], [-4, 4], 'mix' , ...
             [-2.1012, -2.0469; ...
              -2.6718, -2.9028], ...
             [-0.9288, -0.9659; ...
              -2.1180, -2.1475], true}, ...
            {3, 'point_exact', @RefSol.getPoisson2D, [-1, 1], [-1, 1], 'pure', ...
             [-1.1903, -1.0846; ...
              -1.1049, -1.0449], ...
             [ 0.0718,  0.0147; ...
              -0.0761, -0.0036], true}, ...
            {4, 'point_exact', @RefSol.getPoisson2D, [-1, 1], [-1, 1], 'mix' , ...
             [-1.1277, -1.0843; ...
              -1.0899, -1.0449], ...
             [ 0.1880,  0.0164; ...
              -0.0263, -0.0036], true}
                  };
    end

    methods (Test)
        
        function Basis(~)
            % Tests on basis (function) definitions.
            % TODO: implement. 
        end
        
        function Interpolation(~)
            % Tests on interpolation operator.
            % TODO: implement. 
        end
        
        function Assembling(~)
            % Tests on assembling of mass, stiffnes matrix and rhs vector.
            % TODO: implement. 
        end
        
        function Mapping(~)
            % Tests on DOF mapping.
            % TODO: implement. 
        end

        function Convergence(test_case, factory)
            % Convergence test for FE solution of the Laplace problem.
            
            %% Get/Set test case variables.
            
            % Set (fixed) parameter.
            mesh_type = 'cube';
            ref_steps = 1:4;
            order = [1, 2];
            
            % Get (variable) parameter from input.
            TX = struct();
            [factory_num, TX.type, refSolFun, x, y, bnd_variant, ...
             rates_L2_ref, rates_H1_ref, verbosity] = factory{:};
            
            % Expand input variables.
            switch TX.type
                case 'reference'
                    % Ref. solution.
                    TX.ref_sol_u = refSolFun();
                    % Rhs ref. sol.
                    TX.ref_sol.f = @(X, Y) TX.ref_sol_u.L(X, Y) + TX.ref_sol_u.f(X, Y);
                    
                case 'point_exact'
                    % TX coords and strength.
                    TX.coo = [0, 0];
                    TX.val = 1;
                    % Ref. solution.
                    TX.ref_sol_u = refSolFun(TX.coo);
            end
            
            % RX positions (diagonal profile with 17 points).
            RX = [linspace(x(1), x(end), 17).', ...
                  linspace(y(end), y(1), 17).'];
            
            %% Define boundary conditions.
            
            bnd = struct();
            bnd.name = {'xmin', 'xmax', 'ymin', 'ymax'};
            switch bnd_variant
                case 'pure'
                % Pure Dirichlet BC.
                    bnd.type = {'dirichlet'};
                    %                     xmin            xmax            ymin            ymax
                    bnd.val = {{TX.ref_sol_u.f; TX.ref_sol_u.f; TX.ref_sol_u.f; TX.ref_sol_u.f}};
                    bnd.name = {'xmin', 'xmax', 'ymin', 'ymax'};
                case 'mix'
                    % Dirichlet/Neumann BC.
                    bnd.type = {'neumann', 'dirichlet'};
                    %                     xmin            xmax            ymin            ymax                            
                    bnd.val = {{TX.ref_sol_u.J;             [];             []; TX.ref_sol_u.J}, ...
                               {[];             TX.ref_sol_u.f; TX.ref_sol_u.f;             []}};
                    bnd.quad_ord = TX.ref_sol_u.quad_ord;
            end
            
            %% Set up loops.
                            
            % Iterate (if required) over different Lagrange orders.
            [err_L2, err_H1, err_num_DOF] = deal(cell(length(order), 1));
            for cur_order = order

                % Preallocate variables.
                [err_L2{cur_order}, err_H1{cur_order}, err_num_DOF{cur_order}] = ...
                    deal(zeros(length(ref_steps), 1));

                % Iterate over different uniform refinement steps.
                for cur_ref = ref_steps

                    % Print out test case information.
                    if verbosity
                       fprintf('--------------------------------------\n');
                       fprintf(sprintf('Run Test %d for:\n', factory_num));
                       fprintf(sprintf('- "%d" oder Lagrange elements\n', cur_order));
                       fprintf(sprintf('- "%s" basic mesh\n', mesh_type));
                       fprintf(sprintf('- "%d" mesh refinements\n', cur_ref));
                       fprintf(sprintf('- "%s" source\n', TX.type));
                       if length(bnd.type) > 1
                           fprintf('- mixed boundary conditions\n\n');
                       else
                           fprintf(sprintf('- "%s" boundary conditions\n\n', bnd.type{1}));
                       end
                    end

                    %% Init mesh.

                    mesh = Mesh.initMesh(mesh_type, 'bnd', [x, y], ...
                                         'ref', cur_ref, 'verbosity', verbosity);

                    %% Set up Parameter.

                    % Set const. background parameter for grid.
                    param = 1 + zeros(length(mesh.cell2vtx), 1);

                    %% Set up FE structure.

                    % Try to set up interpolation matrix in only one case.
                    if cur_ref == 1
                        cur_RX = RX;
                    else
                        cur_RX = [];
                    end
                    fe = FeL.initFiniteElement(cur_order, mesh, cur_RX, verbosity);

                    %% Set up BC.

                    cur_bnd = FeL.assignBC(bnd, fe, mesh, param);

                    %% Set up FEM linear System.

                    % Set up system matrix.
                    % (for Poisson/Laplace, this only comprises the stiffness matrix)
                    sol.A = FeL.assembleStiff(fe, mesh, param, verbosity);
                    if strcmp(TX.type, 'reference')
                        sol.A = sol.A + FeL.assembleMass(fe, mesh, param, verbosity);
                    end

                    % Set up rhs vector.
                    sol.b = FeL.assembleRHS(fe, mesh, TX, verbosity);

                    % Handle boundary conditions.
                    sol = FeL.treatBC(fe, mesh, sol, cur_bnd, verbosity);

                    if verbosity
                       fprintf('... Linear system and BC set up.\n \n'); 
                    end

                    %% Solve fwd problem.

                    % Get solution at DOF.
                    u = FeL.solveFwd(sol, fe, verbosity);

                    %% Calculate errors.

                    cur_idx = cur_ref - ref_steps(1) + 1;
                    [err_L2{cur_order}(cur_idx), err_H1{cur_order}(cur_idx)] = ...
                        Test.getError(mesh, fe, u, TX.ref_sol_u);
                    err_num_DOF{cur_order}(cur_idx) = fe.sizes.DOF;
                end
            end

            %% Show results.

            % Print  statistics.
            for kk = 1:length(order)
                fprintf(sprintf('Order %d: \n Num DOF  L2 error \t H1 error \n', ...
                    order(kk)));
                for ii = 1:length(ref_steps)
                    fprintf(sprintf('%d \t %d \t %d \n', ...
                        err_num_DOF{kk}(ii), err_L2{kk}(ii), err_H1{kk}(ii)));
                end
            end
            [fit_L2, last_rate_L2] = Test.getConvRates('L2', err_num_DOF, err_L2, true);
            [fit_H1, last_rate_H1] = Test.getConvRates('H1', err_num_DOF, err_H1, true);

            % Plot convergence rates.
            fig_name = sprintf(['" %s " source with " ', ...
                repmat('%s ', 1, length(bnd.type)), '" boundary condition'], ...
                    TX.type, bnd.type{:});
            figure('Name', fig_name, 'NumberTitle', 'off');
            for kk = 1:length(order)
               subplot(2, 1, kk)
                   h_x = logspace(log10(err_num_DOF{kk}(1)), ref_steps(end) + 1, 100);
                   h = fliplr(h_x);
                   loglog(err_num_DOF{kk}, err_L2{kk}/err_L2{kk}(1), 'x-', ...
                       err_num_DOF{kk}, err_H1{kk}/err_H1{kk}(1), 'o-', ...
                       h_x, ones(size(h)), '-b', ...
                       h_x, h/h(1), '-k' ...
                       );
                   x_lim = vertcat(err_num_DOF{:});
                   x_lim = [10^(floor(log10(min(x_lim)))), ...
                       10^(ceil(log10(max(x_lim))))];
                   xlim(x_lim);
                   ylim([1e-3, 1e1]);
                   title(sprintf('u_{FE} vs. u_{ref} w.r.t. DOF (%d. order Lagrange)', ...
                       order(kk)));
                   ylabel('error');
                   xlabel('DOF');
                   legend('||u_{FE} - u_{ref}||_{L2}', '||u_{FE} - u_{ref}||_{H1}', ...
                       'const', 'O(h)', ...
                       'Location', 'EastOutside');
            end
            
            %% Compare to known results.
            
            % Summarize convergence rates.
            rates_L2 = cell2mat([fit_L2, last_rate_L2]);
            rates_H1 = cell2mat([fit_H1, last_rate_H1]);
            
            % Test against known rates.
            test_case.verifyEqual(rates_L2, rates_L2_ref, 'AbsTol', 1e-3);
            test_case.verifyEqual(rates_H1, rates_H1_ref, 'AbsTol', 1e-3); 
        end
    end

end