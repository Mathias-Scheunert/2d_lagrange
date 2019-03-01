function eleVD_fun = getElectrodeAtVertDike(rho1, rho2, x0, h, TX, I)
    % Semi-analytic solution for point source at HS with a vertical dike.
    %
    % SYNTAX
    %   U = getElectrodeAtVertDike(rho1, rho2, x0, h, TX, I)
    %
    % INPUT PARAMETER
    %   rho1 ... Resistivity of surrounding half-space.
    %   rho2 ... Resistivity of vertical dike.
    %   x0   ... Location of dike left boundary.
    %   h    ... Thickness of dike.
    %   TX   ... Vector of source position.
    %   I    ... Electric current.
    %
    % OUTPUT PARAMETER
    %   eleVD_fun ... Struct, containing function handle.
    %                 Handle evaluates Potential at [X, Y = 0]
    %
    % REMARKS
    % Formulas are derived analog to the formulas of Telford et al. 
    % (Applied Geophysics p.683) and Keller & Frischknecht (Electromagnetic
    % Methods in Geophysical Prospecting p.181). 
    % However the formulas are not the same, the formulas in this program 
    % are more general. 
    % The summation is done with the Horner schema. 
    % The survey configuration is a combination of this current-potential
    % configuration.
    %
    % Dike in half-space:          Telford et al, 2nd Ed., 1990 (p. 574)
    % Dike between quarter-spaces: Roy K.K., 2007 (p. 334)
    %
    % COORDINATE SYSTEM
    %
    %                   x = strike direction          
    %                  /            /
    %                 /            /
    %                y = 0        y = h
    %  x,z=0 --------+------------+----------> y
    %         rho1   |    rho2    |    rho1
    %                |            |
    %                v            
    %                z
    %
    % Note the difference between the local coordinates and INPUT:
    %   y -> x for TX
    %              
    % COPYRIGHT
    %   Function was originally developed by T. Hanstein.
    
    %% Check input.
    
    assert(isscalar(rho1), 'rho1 - Scalar denoting resistivity, expected.');
    assert(isscalar(rho2), 'rho2 - Scalar denoting resistivity, expected.');
    assert(isscalar(x0), 'y0 - Scalar denoting dike left boundary, expected.');
    assert(isscalar(h), 'h - Scalar denoting first dike thickness, expected.');
    assert(isscalar(I), 'I - Scalar denoting source current, expected.');
    assert(isvector(TX) && length(TX) == 2, ...
        'TX - Vector [2 x 1], denoting source position, expected.');
    assert(TX(2) == 0, ...
        'Sources only at y = 0 are supported.');
    
    %% Set Parameter.
    
    max = 100; % fix number of summations
    pi2 = 2*pi;
    h2 = 2*h;
    AK = (rho2 - rho1) / (rho2 + rho1);
    AK2 = AK*AK;
    CC = 1 + AK;
    DD = 1 - AK2;
    AA = -AK * DD;
    BB = -AK * CC;

    %% Define function.
    
    % (Note the change of coordinates.)
    eleVD_fun.f = @(X, Y) evaluateU(0, (TX(1) - x0), 0, (X - x0), Y);
    
    %% Define gradient.
    
    % Skip derivation.

    %% Set required quadrature order.
    
    eleVD_fun.quad_ord = 4;

    %% Evaluate cases for different positions of TX and RX.
    
    % XC ... x-coord. of TX, parallel to strike direction
    % YC ... y-coord. of TX, perpendicular to strike direct
    % XP ... x-coord. of RX, parallel to strike direction
    % YP ... y-coord. of RX, perpendicular to strike direct
    
    function U = evaluateU(XC, YC, XP, YP, ZP)
        assert(ZP == 0, ...
            'Evaluation only at y = 0 supported.');
        X = XP - XC;
        X2 = X*X;
        if (XC == XP) && (YC == YP)
            % TX and RX coincide.
            U = 0;
            return
        elseif YC <= 0
            % TX left of dike
            S = abs(YC);
            S2 = S + S;
            if YP <= 0
                % RX left of dike
                A = YP - YC;
                Y = (max+1)*h2 + S2 - A;
                SUM = 1/sqrt(X2 + Y*Y);
                for mm = max:-1:1
                    Y = mm*h2 + S2 - A;
                    SUM = 1/sqrt(X2 + Y*Y) + AK2*SUM;
                end
                V = 1/sqrt(X2+A*A) + AK/sqrt(X2 + (S2-A).^2) + AA*SUM;
                U = (I*rho1/pi2) * V;
                return
            elseif YP <= h
                % RX within dike
                A = abs(YP) + abs(YC);
                Y = (max+1)*h2 + S2 - A;
                SUM1 = 1/sqrt(X2 + Y*Y);
                Y = max*h2 + A;
                SUM2 = 1/sqrt(X2 + Y*Y);
                for mm = max:-1:1
                    Y = mm*h2 + S2 - A;
                    SUM1 = 1/sqrt(X2 + Y*Y) + AK2*SUM1;
                    Y = (mm-1)*h2 + A;
                    SUM2 = 1/sqrt(X2 + Y*Y) + AK2*SUM2;
                end
                V = BB*SUM1 + CC*SUM2;
                U = (I*rho1/pi2) * V;
                return
            else
                % RX right of dike
                A = abs(YP) + abs(YC);
                Y = max*h2 + A;
                SUM = 1/sqrt(X2 + Y*Y);
                for mm = max:-1:1
                    Y = (mm-1)*h2 + A;
                    SUM = 1/sqrt(X2 + Y*Y) + AK2*SUM;
                end
                V = DD*SUM;
                U = (I*rho1/pi2) * V;
                return
            end
        elseif (YC > 0) && (YC <= h)
            % TX within dike
            S = YC;
            S2 = S + S;
            if YP <= 0
                % RX left of dike
                A = abs(YP) + YC;
                Y = max*h2 + A;
                SUM1 = 1/sqrt(X2 + Y*Y);
                Y = (max+1)*h2 - S2 + A;
                SUM2 = 1/sqrt(X2 + Y*Y);
                for mm = max:-1:1
                    Y = (mm-1)*h2 + A;
                    SUM1 = 1/sqrt(X2 + Y*Y) + AK2*SUM1;
                    Y = mm*h2 - S2 + A;
                    SUM2 = 1/sqrt(X2 + Y*Y) + AK2*SUM2;
                end
                V = (1 + AK)*( SUM1 - AK*SUM2);
                U = (I*rho1/pi2) * V;
                return
            elseif YP <= h
                % RX within dike
                A = YP - YC;
                Y = (max+1)*h2 - A;
                SUM1 = 1/sqrt(X2 + Y*Y);
                Y = (max+1)*h2 + A;
                SUM2 = 1/sqrt(X2 + Y*Y);
                Y = max*h2 + S2 + A;
                SUM3 = 1/sqrt(X2 + Y*Y);
                Y = (max+1)*h2 - S2 - A;
                SUM4 = 1/sqrt(X2 + Y*Y);
                for mm = max:-1:1
                    Y = mm*h2 - A;
                    SUM1 = 1/sqrt(X2 + Y*Y) + AK2*SUM1;
                    Y = mm*h2 + A;
                    SUM2 = 1/sqrt(X2 + Y*Y) + AK2*SUM2;
                    Y = (mm-1)*h2 + S2 + A;
                    SUM3 = 1/sqrt(X2 + Y*Y) + AK2*SUM3;
                    Y = mm*h2 - S2 - A;
                    SUM4 = 1/sqrt(X2 + Y*Y) + AK2*SUM4;
                end
                V = 1/sqrt(X2+A*A) + AK2*(SUM1+SUM2) - AK*(SUM3 + SUM4);
                U = (I*rho2/pi2) * V;
                return
            else
                % RX right of dike
                A = YP - YC;
                Y = max*h2 + A;
                SUM1 = 1/sqrt(X2 + Y*Y);
                Y = max*h2 + S2 + A;
                SUM2 = 1/sqrt(X2 + Y*Y);
                for mm = max:-1:1
                    Y = (mm-1)*h2 + A;
                    SUM1 = 1/sqrt(X2 + Y*Y) + AK2*SUM1;
                    Y = (mm-1)*h2 + S2 + A;
                    SUM2 = 1/sqrt(X2 + Y*Y) + AK2*SUM2;
                end
                V = (1 + AK)*(SUM1 - AK*SUM2);
                U = (I*rho1/pi2) * V;
                return
            end
        elseif YC > h
            % TX right of dike
            S = YC - h;
            S2 = S + S;
            if YP <= 0
                % RX left of dike
                A = abs(YP) + YC;
                Y = max*h2 + A;
                SUM = 1/sqrt(X2 + Y*Y);
                for mm = max:-1:1
                    Y = (mm-1)*h2 + A;
                    SUM = 1/sqrt(X2 + Y*Y) + AK2*SUM;
                end
                V = DD*SUM;
                U = (I*rho1/pi2) * V;
                return
            elseif YP <= h
                % RX within dike
                A = abs(YP - YC);
                Y = max*h2 + A;
                SUM1 = 1./sqrt(X2 + Y*Y);
                Y = (max + 1)*h2 + S2 - A;
                SUM2 = 1/sqrt(X2 + Y*Y);
                for mm = max:-1:1
                    Y = (mm-1)*h2 + A;
                    SUM1 = 1/sqrt(X2 + Y*Y) + AK2*SUM1;
                    Y = mm*h2 + S2 - A;
                    SUM2 = 1/sqrt(X2 + Y*Y) + AK2*SUM2;
                end
                V = CC*(SUM1 - AK*SUM2);
                U = (I*rho1/pi2) * V;
                return
            else
                % RX right of dike
                A = YP - YC;
                Y = (max+1)*h2 + S2 + A;
                SUM = 1/sqrt(X2 + Y*Y);
                for mm = max:-1:1
                    Y = mm*h2 + S2 + A;
                    SUM = 1/sqrt(X2 + Y*Y) + AK2*SUM;
                end
                V = 1/sqrt(X2 + A*A) + AK/sqrt(X2 + (S2 + A).^2) + AA*SUM;
                U = (I*rho1/pi2) * V;
                return
            end
        end
    end
end