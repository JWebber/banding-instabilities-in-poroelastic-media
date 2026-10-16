function [phi, a, q] = base_state(phiC, alpha, M, P)
    % Solve for \Phi(0)
    basal_fun = @(x)((x-1)./(1-x/phiC) + 4*M*(x.^(1/3)-1)-P);
    ops = optimset('display', 'none');
    basal_phi = fsolve(basal_fun, 1, ops);

    % DE is (diff*\Phi')' = 0, 
    diffusivity = @(y)(y(1)^(-alpha).*((1-1/phiC)./(1-y(1)/phiC).^2 + 4*M*y(1).^(-2/3)/3));
    eqn = @(t, y)[y(2)./diffusivity(y); 0];
    bvs = @(ya, yb)([ya(1)-basal_phi; yb(1)-1]);
    bv_ops = bvpset('RelTol', 1e-8);
    system = bvpinit(linspace(0, 1, 250), @(z)[basal_phi + (1-basal_phi)*z; 0]);
    sol = bvp5c(eqn, bvs, system, bv_ops);
    
    raw_phi = @(z)interp1(sol.x, sol.y(1, :), z);
    a = 1/integral(raw_phi, 0, 1);

    phi = @(z)interp1(sol.x, sol.y(1, :), z/a);
    q = (sol.y(2, 1))/a;
end

