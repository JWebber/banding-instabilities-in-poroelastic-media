% Plot the dev strain in the growing instability to illustrate the planar mechanism
M = 0.1;
alpha = 2;

guess = @(z)[1;0;z;1;0;1;1];
guess = bvpinit(linspace(0,1,100), guess, 1);

[s, guess, thickness] = solve_zero_gr(20, M, alpha, guess);
phi1 = -guess.y(1, :).^(2/3).*guess.y(4, :);
p = @(z)interp1(thickness*guess.x, guess.y(1, :), z);
eta1p = @(z)interp1(thickness*guess.x, guess.y(4, :), z);
newPhi = @(epsi, u)(interp1(guess.x, guess.y(1, :) + epsi*phi1, u));
newthk = @(epsi)(1/integral(@(u)newPhi(epsi, u), 0, 1));

mat_write = zeros(100, 20);

i = 1;
for epsi = linspace(0, 0.5, 10)
    epszz = @(z, epsi)(2*(1-p(z).^(1/3))+(2/3)*epsi*eta1p(z));
    plot(newthk(epsi)*linspace(0, 1, 100), epszz(linspace(0, thickness, 100), epsi));
    mat_write(:, i) = newthk(epsi)*linspace(0, 1, 100)';
    mat_write(:, i+1) = epszz(linspace(0, thickness, 100), epsi)';

    i = i+2;
    hold on;
end

writematrix(mat_write, 'output/mechanism.csv');