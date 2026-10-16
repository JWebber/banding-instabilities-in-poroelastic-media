function [s, sol, thickness] = solve_zero_gr(pressureHead, M, alpha, solinit)
    diffusivity = @(p)(p.^(-alpha).*(1 + 4*M*p.^(-2/3)/3));
    diffusivityP = @(p)(-alpha*p.^(-alpha-1) - (2/3+alpha)*(4*M/3)*p.^(-(alpha+5/3)));
    ops = optimset('display', 'none');

    % SETUP FOR EQUATIONS; 'P' denotes 'prime' so pPEqn is the equation for
    % (p')' (or \Phi_0''); equations are stated for \Phi_0, Z, \eta_1 in
    % that order
    % ---

    % \Phi
    pPEqn = @(p, pP)(-(diffusivityP(p)./diffusivity(p)).*pP.^2);

    % \eta
    prefactor = @(p)(p.^(2/3) + 4*M/3);
    Q = @(p, pP)(diffusivity(p).*pP);
    bkt1 = @(s, p, pP)(s.*p.^(alpha-1/3) - (alpha/3)*Q(p, pP).*p.^(alpha-4/3).*pP - (4/9)*p.^(-4/3).*pP.^2 - (2/3)*p.^(-1/3).*pPEqn(p, pP));
    bkt2 = @(p, pP)(alpha*Q(p, pP).*p.^(alpha-1/3) - (4*M/3)*pP./p - (7/3)*p.^(-1/3).*pP);
    etaPPEqn = @(e, eP, ePP, s, p, pP)((bkt1(s, p, pP).*eP + bkt2(p, pP).*ePP - s.*p.^(alpha-4/3).*pP.*e)./prefactor(p));
    
    % ---
    % SOLUTION VECTOR IS [\Phi_0; \Phi_0'; \eta_1; \eta_1'; \eta_1''; N; I]
    % ---

    % Build the system as a vector, with the auxiliary normalisation
    % function N = \int{\eta^2}
    eqSystem = @(z, y, s)[y(2); pPEqn(y(1), y(2)); y(4); y(5); etaPPEqn(y(3), y(4), y(5), s, y(1), y(2)) ; y(3).^2 ; y(1)];

    % SETUP FOR BCs: seventh-order system w/ one parameter requires eight BCs
    % ---

    % \Phi BCs (2)
    basalPhiFun = @(z)(z-1 + 4*M*(z.^(1/3)-1) - pressureHead);
    if(pressureHead ~= 0 && basalPhiFun(1)*basalPhiFun(1+pressureHead)<0)
        basalPhi = fzero(basalPhiFun, [1, 1+pressureHead], ops);
    elseif(pressureHead ~= 0)
        basalPhi = fsolve(basalPhiFun, 1+pressureHead, ops);
    else
        basalPhi = 1;
    end

    pBCs = @(pa, pb)[pb-1; pa - basalPhi]; % \Phi = 1 on top, plus base val

    % Pert BCs (3)
    pertBCs = @(eBase, ePBase, eTop, ePTop, pPTop)[eBase; ePBase; Q(1, pPTop)*eTop - (1+4*M/3)*ePTop];

    % Auxiliary BCs (3)
    auxBCs = @(Na, Nb, Ia)[Na; Nb-1; Ia];

    bcSystem = @(ya, yb, s)[pBCs(ya(1), yb(1)); pertBCs(ya(3), ya(4), yb(3), yb(4), yb(2)); auxBCs(ya(6), yb(6), ya(7))];
    

    % SOLVER
    % ---

    options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8, 'Nmax', 40000);
    sol = bvp5c(eqSystem, bcSystem, solinit, options);
    s = sol.parameters(1);

    thickness = 1/sol.y(7, end);

end