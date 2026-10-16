% SWEEP DIFFERENT VALUES OF P
M = 1;
alpha = 2/3;
phiC = 10;

Y = zeros(100, 100);
P = logspace(-1, 2, 100);
X = repmat(P, [100, 1]);

Z = zeros(100, 100);

for i = 1 : 100
    [p, a, ~] = base_state(phiC, alpha, M, P(i));
    Y(:, i) = linspace(0, a, 100)';
    Z(:, i) = p(linspace(0, a, 100)');
end

contourf(X, Y, Z, 250, 'edgecolor', 'none');
set(gca, 'xscale', 'log');
colorcet('L20');
clim([1 10]);
fprintf(strcat(num2str(max(Y(:))), '\n')); % output the thickness
high_quality_png('output/varyP.png');

% SWEEP DIFFERENT VALUES OF M
P = 50;
M = logspace(-2, 2, 100);
X = repmat(M, [100, 1]);

for i = 1 : 100
    [p, a, ~] = base_state(phiC, alpha, M(i), P);
    Y(:, i) = linspace(0, a, 100)';
    Z(:, i) = p(linspace(0, a, 100)');
end

contourf(X, Y, Z, 250, 'edgecolor', 'none');
set(gca, 'xscale', 'log');
colorcet('L20');
clim([1 10]);
fprintf(strcat(num2str(max(Y(:))), '\n')); % output the thickness
high_quality_png('output/varyM.png');

% SWEEP DIFFERENT VALUES OF ALPHA
M = 1;
alpha = linspace(0, 5, 100);
X = repmat(alpha, [100, 1]);

for i = 1 : 100
    [p, a, ~] = base_state(phiC, alpha(i), M, P);
    Y(:, i) = linspace(0, a, 100)';
    Z(:, i) = p(linspace(0, a, 100)');
end

contourf(X, Y, Z, 250, 'edgecolor', 'none');
set(gca, 'xscale', 'linear');
colorcet('L20');
clim([1 10]);
fprintf(strcat(num2str(max(Y(:))), '\n')); % output the thickness
high_quality_png('output/varyA.png');

% SWEEP DIFFERENT VALUES OF PHIC
alpha = 2/3;
phiC = linspace(1.5, 10, 100);
X = repmat(phiC, [100, 1]);

for i = 1 : 100
    [p, a, ~] = base_state(phiC(i), alpha, M, P);
    Y(:, i) = linspace(0, a, 100)';
    Z(:, i) = p(linspace(0, a, 100)');
end

contourf(X, Y, Z, 250, 'edgecolor', 'none');
colorcet('L20');
set(gca, 'xscale', 'linear');
clim([1 10]);
fprintf(strcat(num2str(max(Y(:))), '\n')); % output the thickness
high_quality_png('output/varyPC.png');