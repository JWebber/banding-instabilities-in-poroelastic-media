N = 100;
alpha = 2;

gamma_range = [0 logspace(-2, log10(4), N-1)];

matrix = zeros(N, 11);
maxmatrix = zeros(10, 2);

matrix(:, 1) = gamma_range';

% i = 2;
% 
% for M = linspace(0.01, 1, 10)
%     ratesI = get_growth_rates(15, M, 2, gamma_range);
%     plot(gamma_range, ratesI);
%     matrix(:, i) = ratesI';
%     i = i+1;
%     hold on;
%     [val, mp] = max(ratesI);
%     plot(gamma_range(mp), val, 'r.');
%     drawnow;
%     maxmatrix(i-1, :) = [gamma_range(mp), val];
% end
% 
% writematrix(matrix, 'output/vary_m.csv');
% writematrix(maxmatrix, 'output/vary_m_max.csv');

i = 2;

for alpha = linspace(1.5, 3, 10)
    ratesI = get_growth_rates(15, 0.75, alpha, gamma_range);
    plot(gamma_range, ratesI);
    hold on;
    matrix(:, i) = ratesI';
    i = i+1;
    [val, mp] = max(ratesI);
    plot(gamma_range(mp), val, 'r.');
    drawnow;
    maxmatrix(i-1, :) = [gamma_range(mp), val];
end

writematrix(matrix, 'output/vary_a.csv');
writematrix(maxmatrix, 'output/vary_a_max.csv');

function rates = get_growth_rates(P, M, alpha, gammas)
    % Initialise the guess to start
    guess = @(z)[1;0;z;1;0;1;1];
    guess = bvpinit(linspace(0,1,100), guess, 1);

    % Jump forward to the value of P we are considering (retain the right
    % branch!)
    for pJump = linspace(min(P/5, 1e-1), P, 10)
        [s0, guess, t0] = solve_zero_gr(pJump, M, alpha, guess);
    end

    sz = max(size(gammas, 1), size(gammas, 2));
    rates = NaN*zeros(1, sz);

    finite_gamma_guess = transform_guess(guess);

    for i = 1 : sz
        % Two separate cases; if we need the zero-gamma solution start there
        if(gammas(i) == 0)
            rates(i) = s0/t0^2;
        else
            % Don't forget to convert gamma_real to gamma_normalised here;
            % hence multiply by t0
            [s, finite_gamma_guess, ~] = solve_bvp(P, M, alpha, gammas(i)*t0, finite_gamma_guess);
            rates(i) = s/t0^2;
        end
    end
end