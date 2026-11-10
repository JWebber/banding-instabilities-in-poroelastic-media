guess = @(z)[2-z;-1;z.^2;2*z;2;1;1];
guess = bvpinit(linspace(0,1,100), guess, 1);
alpcrit = 2;
[~, guess, ~] = solve_zero_gr(10, 0.01, alpcrit, guess);

for M = [0.01 0.1 0.25 0.5 0.75 1 1.25 1.5 1.75 2]%[0.01 0.1 0.25 0.5 0.75 1 2.5 5 10]
    matrix = [logspace(log10(1.5), log10(10), 40)' zeros(40, 1); logspace(log10(10), log10(1000), 40)' zeros(40, 1)];
    
    guesspivot = guess;
    alppivot = alpcrit;

    j = 40;
    
    for P = fliplr(logspace(log10(1.5), log10(10), 40))
        alpcrit = fsolve(@(x)pureGr(x, P, M, guess), alpcrit);
        [~, guess, ~] = solve_zero_gr(P, M, alpcrit, guess);
    
        plot(P, alpcrit, 'r.');
        hold on;
        ylim([0 5]);
        drawnow;

        matrix(j, 2) = alpcrit;
        j = j-1;

        if(alpcrit > 5)
            break;
        end
    end
    
    alpcrit = alppivot;
    guess = guesspivot;
    j = 41;
    
    for P = logspace(log10(10), log10(1000), 40)
        alpcrit = fsolve(@(x)pureGr(x, P, M, guess), alpcrit);
        [~, guess, ~] = solve_zero_gr(P, M, alpcrit, guess);
    
        plot(P, alpcrit, 'r.');
        hold on;
        ylim([0 5]);
        drawnow;

        matrix(j, 2) = alpcrit;
        j = j+1;
    end

    writematrix(matrix, strcat('output/crit_a_m_', num2str(M), '.csv'));

    guess = guesspivot;
    alpcrit = alppivot;

end

function s = pureGr(alpha, P, M, guess)
    [s, ~, ~] = solve_zero_gr(P, M, alpha, guess);
end