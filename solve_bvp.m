function [s, sol, thickness] = solve_bvp(pressureHead, M, alpha, gamma, solinit, varargin)
    diffusivity = @(p)(p.^(-alpha).*(1 + 4*M*p.^(-2/3)/3));
    diffusivityP = @(p)(-alpha*p.^(-alpha-1) - (2/3+alpha)*(4*M/3)*p.^(-(alpha+5/3)));
    ops = optimset('display', 'none');

    % SETUP FOR EQUATIONS; 'P' denotes 'prime' so pPEqn is the equation for
    % (p')' (or \Phi_0''); equations are stated for \Phi_0, Z, \eta_1 and \Psi in
    % that order
    % ---

    % \Phi
    pPEqn = @(p, pP)(-(diffusivityP(p)./diffusivity(p)).*pP.^2);
    
    % Z
    ZPEqn = @(Z)(gamma^2.*Z);

    % \eta_1
    etaPEqn = @(e, ZP, PP)(gamma^2.*e + ZP/M - PP);

    % \Psi
    prefactor = @(p)(p.^(2/3) + 4*M/3);
    Q = @(p, pP)(diffusivity(p).*pP);
    bkt1 = @(s, p, pP)(s.*p.^(alpha-1/3) - (alpha/3)*Q(p, pP).*p.^(alpha-4/3).*pP - (4/9)*p.^(-4/3).*pP.^2 - (2/3)*p.^(-1/3).*pPEqn(p, pP));
    bkt2 = @(p, pP)(alpha*Q(p, pP).*p.^(alpha-1/3) - (4*M/3)*pP./p - (7/3)*p.^(-1/3).*pP);
    PPEqn = @(P, PP, e, s, p, pP)(gamma^2*P + (bkt1(s, p, pP).*P + bkt2(p, pP).*PP + s.*p.^(alpha-4/3).*pP.*e)./prefactor(p));
    
    % ---
    % SOLUTION VECTOR IS [\Phi_0; \Phi_0'; Z; Z'; \eta_1; \eta_1'; \Psi;
    % \Psi'; N; I]
    % ---

    N = @(z, y)(y(3).^2 + y(7).^2 + y(5).^2 + ((y(7) + y(6))/gamma).^2);

    if(nargin == 6) % Use the linearised case
        etaBase = @(z)interp1(solinit.x, solinit.y(5, :), z);
        ZBase = @(z)interp1(solinit.x, solinit.y(3, :), z);
        PBase = @(z)interp1(solinit.x, solinit.y(7, :), z);
        xiBase = @(z)interp1(solinit.x, (solinit.y(7, :) + solinit.y(6, :))/gamma, z);

        N = @(z, y)(y(3).*ZBase(z) + y(7).*PBase(z) + y(5).*etaBase(z) + xiBase(z).*(y(7, :) + y(6, :))/gamma);
    end

    % Build the system as a vector, with the auxiliary normalisation
    eqSystem = @(z, y, s)[y(2); pPEqn(y(1), y(2)); y(4); ZPEqn(y(3)); y(6); etaPEqn(y(5), y(4), y(8)); y(8); PPEqn(y(7), y(8), y(5), s, y(1), y(2)); N(z, y); y(1)];

    % SETUP FOR BCs: tenth-order system w/ one parameter requires eleven BCs
    % ---

    % \Phi BCs (2)
    basalPhiFun = @(z)(z-1 + 4*M*(z.^(1/3)-1) - pressureHead);
    if(pressureHead ~= 0)
        basalPhi = fzero(basalPhiFun, [1, 1+pressureHead], ops);
    else
        basalPhi = 1;
    end

    pBCs = @(pa, pb)[pb-1; pa - basalPhi]; % \Phi = 1 on top, plus base val

    % Base BCs (3)
    baseBCs = @(Z, e, eP, P)[e; P + eP; Z - (basalPhi^(2/3) + 4*M/3)*P];

    % Top BCs (3)
    topBCs = @(Z, ZP, e, eP, P, pP)[2*M*(eP+P)-Z; ZP + 2*M*gamma^2*e; Z - (1 + 4*M/3)*P - e*Q(1, pP)];

    % Auxiliary BCs (3)
    auxBCs = @(Na, Nb, Ia)[Na; Nb-1; Ia];

    bcSystem = @(ya, yb, s)[pBCs(ya(1), yb(1)); baseBCs(ya(3), ya(5), ya(6), ya(7)); topBCs(yb(3), yb(4), yb(5), yb(6), yb(7), yb(2)); auxBCs(ya(9), yb(9), ya(10))];
    

    % SOLVER
    % ---

    options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8, 'Nmax', 40000);
    sol = bvp5c(eqSystem, bcSystem, solinit, options);
    s = sol.parameters(1);

    thickness = 1/sol.y(10, end);
end