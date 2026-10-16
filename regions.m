% Key parameters
alpha = 2;

% Resolution of M
Nm = 50;
Ms = logspace(log10(0.05), 1, Nm);

region_starts = NaN*zeros(1, Nm);
region_ends = NaN*zeros(1, Nm);
region_crossovers = NaN*zeros(1, Nm);

% Initialise the guess to start
guess = @(z)[1;0;z;1;0;1;1];
guess_end_0 = bvpinit(linspace(0,1,100), guess, 1);

% Where do we start looking?
guess_start = 10.5531; %37; %10.5531;
guess_end = 10.5903; %38; %10.5903; 
guess_crossover_v = 81;
ops = optimset('tolx', 5e-3);

t = tiledlayout(1, 2, 'tilespacing', 'loose', 'padding', 'tight');
nexttile;
plot([guess_start guess_end], [0 0]);
a1 = gca;
nexttile(t);
plot([80 81], [0 0]);
a2 = gca;

for i = 1 : Nm
    % Start by finding the upper bound where s(0)=0

    % Initialise ahead of the first step
    for P = linspace(1e-1, guess_end, 5)
        [s, guess_end_0, thickness] = solve_zero_gr(P, Ms(i), alpha, guess_end_0);
    end

    guess_end = fzero(@(z)upperBoundFun(z, Ms(i), alpha, guess_end_0), [0.95*guess_end, 1.05*guess_end], ops);

    
    plot(a1, guess_end, Ms(i), 'r.');
    hold(a1, "on");
    drawnow;

    % Now find the lower bound where s < 0 for all gamma
    guess_start = fzero(@(z)lowerBoundFun(z, Ms(i), alpha, guess_end_0), [0.95*guess_start, 1.05*guess_start], ops);

    % Now find the AB crossover region
    guess_crossover = guess_end_0;
    % Initialise ahead of the first step
    for P = linspace(1e-1, guess_crossover_v, 10)
        [~, guess_crossover, ~] = solve_zero_gr(P, Ms(i), alpha, guess_crossover);
    end
    crossover = fzero(@(z)ABBoundaryFun(z, Ms(i), alpha, guess_crossover), [0.9*guess_crossover_v, 1.1*guess_crossover_v], ops);
    guess_crossover_v = crossover;

    plot(a2, crossover, Ms(i), 'g.');
    hold(a2, "on");
    drawnow;

    region_starts(i) = guess_start;
    region_ends(i) = guess_end;
    region_crossovers(i) = crossover;

    plot(a1, guess_start, Ms(i), 'b.');
    drawnow;

    if (i < Nm && i > 1) % Extrapolate to get the next guesses
        deltaS = guess_start - region_starts(i-1);
        deltaE = guess_end - region_ends(i-1);

        guess_start = guess_start + (Ms(i+1)-Ms(i))*deltaS/(Ms(i)-Ms(i-1));
        guess_end = guess_end + (Ms(i+1)-Ms(i))*deltaE/(Ms(i)-Ms(i-1));
    end

    plot(a1, region_starts, Ms);
    plot(a1, region_ends, Ms);
    plot(a2, region_crossovers, Ms);
    drawnow;
    
    writematrix([Ms', region_starts', region_ends', region_crossovers'], "output/alpha_2_regions.csv")
end

function s = upperBoundFun(P, M, alpha, guess0)
    [s, ~, thickness] = solve_zero_gr(P, M, alpha, guess0);
    s = s / thickness^2;
end

function s = lowerBoundFun(P, M, alpha, guess0)
    [s, g0, thickness] = solve_zero_gr(P, M, alpha, guess0);

    [sN, guess, tN] = solve_bvp(P, M, alpha, 1e-3, transform_guess(g0));
    
    if(sN/tN^2 < s/thickness^2)
        s = s/thickness^2;
        return;
    else
        ops = optimset('tolx', 5e-3);
        sp = fminbnd(@(z)(-1*maxFun(z, P, M, alpha, guess)), 1e-3, 1, ops);
        s = maxFun(sp, P, M, alpha, guess);
    end
end

function s = maxFun(gamma, P, M, alpha, guess)
    [s, ~, t] = solve_bvp(P, M, alpha, gamma, guess);
    s = s/t^2;
end

function v = ABBoundaryFun(P, M, alpha, guess0)
    % Initialise ahead of the first step
    [s, guess0, thickness] = solve_zero_gr(P, M, alpha, guess0);

    [sN, ~, tN] = solve_bvp(P, M, alpha, 1e-2, transform_guess(guess0));

    v = (sN/tN^2) - (s/thickness^2);
end