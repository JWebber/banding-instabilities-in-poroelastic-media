alpha = 2;
Ms = [0.1 0.5 1];
Np = 400;
gams = NaN*zeros(3, Np);
vals = gams;
Ps = linspace(10, 100, Np);

% Take snapshots at different stiffnesses
for i = 1 : 4
    guess0 = @(z)[1;0;z;1;0;1;1];
    guess0 = bvpinit(linspace(0,1,100), guess0, 1);

    gamStar = 1e-6;

    for P = linspace(0.01, min(Ps), 5)
        [s, guess0, t] = solve_zero_gr(P, Ms(i), alpha, guess0);
    end

    [pos, val] = fmincon(@(z)(-maxFun(z, min(Ps), Ms(i), alpha, transform_guess(guess0))), pos, [], [], [], [], 0, 2.5);

    for j = 1 : Np
        [s, guess0, t] = solve_zero_gr(Ps(j), Ms(i), alpha, guess0);

        [pos, val] = fmincon(@(z)(-maxFun(z, Ps(j), Ms(i), alpha, transform_guess(guess0))), pos, [], [], [], [], 0.85*lb, 1.15*pos);

        if(val < 0)
            gams(i, j) = pos/t;
            vals(i, j) = val/t^2;

            if(gams(i, j) < 1/2)
                break;
            end
        end

        plot(Ps, 1./gams');
        ylim([0 1.5]);
        drawnow;
    end
end

function s = maxFun(gamma, P, M, alpha, guess)
    [s, ~, t] = solve_bvp(P, M, alpha, gamma, guess);
    s = s/t^2;
end