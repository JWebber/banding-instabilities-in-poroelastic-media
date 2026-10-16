N = 20;
M = logspace(-2, log10(5), N);
P = logspace(log10(1/2), log10(50), N);
[X, Y] = meshgrid(P, M);

for alpha = [1 2 3]
    growth_rates = NaN*zeros(N, N);
    guess = @(z)[1;0;z;1;0;1;1];
    guess = bvpinit(linspace(0,1,100), guess, 1);
    
    indices = 1 : N;
    
    for i = 1 : N
        for j = indices
            [s, guess, thickness] = solve_zero_gr(P(j), M(i), alpha, guess);
            growth_rates(i, j) = s/thickness^2;
        end
    
        indices = fliplr(indices);
        contourf(X, Y, growth_rates, 400, 'EdgeColor','none');
        set(gca, 'xscale', 'log', 'yscale', 'log');
        clim([-50 50]);
        colorcet('D01A');
        drawnow;
    end
    
    hold on;
    contour(X, Y, growth_rates, [0 0], 'LineWidth', 3, 'edgecolor', [34 139 34]/255);
    
    axis off;
    high_quality_png(strcat('output/zerogamma_alpha_', num2str(alpha), '.png'));
end