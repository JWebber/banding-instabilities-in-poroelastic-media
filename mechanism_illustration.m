% Plot the growing instability to illustrate the planar mechanism
M = 0.1;
alpha = 2;

guess = @(z)[1;0;z;1;0;1;1];
guess = bvpinit(linspace(0,1,100), guess, 1);

[s, guess, thickness] = solve_zero_gr(20, M, alpha, guess);
phi1 = -guess.y(1, :).^(2/3).*guess.y(4, :);

newPhi = @(epsi, u)(interp1(guess.x, guess.y(1, :) + epsi*phi1, u));
newthk = @(epsi)(1/integral(@(u)newPhi(epsi, u), 0, 1));

permeabilityData = NaN*zeros(size(guess.x, 2), 6);
i = 1;

for epsi = [0 0.1 0.3]
    figure;

    X = repmat(linspace(0, 1, 2), [size(guess.x, 2) 1]);
    Y = repmat(newthk(epsi)*guess.x', [1 2]);
    Z = repmat(newPhi(epsi, guess.x)', [1 2]);

    contourf(X, Y, log10(Z-1), 100, 'EdgeColor', 'none');

    ylim([0 1]);
    clim([-1 1]);
    colorcet('L20');
    axis off;

    high_quality_png(strcat('output/planarmech_', num2str(epsi), '.png'));
    
    hold on;

    permeabilityData(:, i) = newthk(epsi)*guess.x';
    permeabilityData(:, i+1) = (newPhi(epsi, guess.x).^(-2))';

    i = i+2;
end

writematrix(permeabilityData, 'output/fig6.csv')