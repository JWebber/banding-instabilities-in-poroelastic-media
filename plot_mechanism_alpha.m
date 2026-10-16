% Vary permeability factor alpha
for alpha = [2 3 4]
    clf;
    guess0 = @(z)[1;0;z;1;0;1;1];
    guess0 = bvpinit(linspace(0,1,100), guess0, 1);

    for P = linspace(0, 12, 10)
        [s, guess0, t] = solve_zero_gr(P, 0.25, alpha, guess0);
    end
    guess = transform_guess(guess0);
    pos = 0.1;
    if(alpha ~= 5)
        [pos, val] = fmincon(@(z)(-maxFun(z, 12, 0.25, alpha, guess)), pos, [], [], [], [], 0, 1);
    
        for gam = linspace(0.05, pos, 10)
            [s, guess, t] = solve_bvp(12, 0.25, alpha, gam, guess);
        end
    else
        guess = transform_guess(guess0);
        pos = 0;
    end

    
    output = interp1(guess.x, guess.y', linspace(0, 1, 8))';

    fprintf(strcat('gam=', num2str(pos/t), ' and s = ', num2str(s/t^2), '\n'));

    phi1 = output(1, :).^(2/3).*output(7, :);
    dim = size(guess.x, 2);

   

    X = repmat(linspace(0, 10, 200), [8, 1]);
    Y = repmat(linspace(0,1,8)', [1 200]);
    Z = NaN*X;
    UX = Z;
    UZ = Z;

    Q = (1 + 2/3)*output(2, end);
    
    p1z = output(4, :) - (2/3 + output(1, :).^(2/3)).*output(8, :) - (2/3)*output(1, :).^(-1/3).*output(2, :).*output(7, :);
    p1 = output(3, :) - (2/3 + output(1, :).^(2/3)).*output(7, :);
    u1z = -output(1, :).^(-alpha).*p1z - alpha*Q*phi1./output(1, :);
    u1x = -output(1, :).^(-alpha).*p1;
    for i = 1 : 200
        Z(:, i) = phi1'*cos(pos*X(1, i)/t);
        UX(:, i) = u1x'.*sin(pos*X(1, i)/t);
        UZ(:, i) = u1z'.*cos(pos*X(1, i)/t);
    end

    % Rescale the arrows so even tiny ones show up with effectively zero
    % length
    mags = UX.^2 + UZ.^2;
    thetas = atan2(UZ, UX);
    smallarrows = (mags < 0.5);

    UX(smallarrows) = 0.5*cos(thetas(smallarrows));
    UZ(smallarrows) = 0.5*sin(thetas(smallarrows));

    contourf(X, Y, Z, 200, 'EdgeColor', 'none');
    clim([-3 3]);
    colorcet('L20');
    axis off;
    high_quality_png(strcat('output/instab_alp_', num2str(alpha), '.png'));
    hold on;
    quiver(X(:, [1:20:200, 200]), Y(:, [1:20:200, 200]), UX(:, [1:20:200, 200]), UZ(:, [1:20:200, 200]));

    ylim([0 1]);

    writematrix([reshape(X(:, [1:20:200, 200]), [], 1), reshape(Y(:, [1:20:200, 200]), [], 1), reshape(UX(:, [1:20:200, 200]), [], 1), reshape(UZ(:, [1:20:200, 200]), [], 1)], strcat('output/quiver_alp_', num2str(alpha), '.csv'))
    writematrix([reshape(X(end, [1:20:200, 200]), [], 1), reshape(Y(end, [1:20:200, 200]), [], 1), reshape(UX(end, [1:20:200, 200]), [], 1), reshape(UZ(end, [1:20:200, 200]), [], 1)], strcat('output/quiver_surf_alp_', num2str(alpha), '.csv'));
end

function s = maxFun(gamma, P, M, alpha, guess)
    [s, ~, t] = solve_bvp(P, M, alpha, gamma, guess);
    s = s/t^2;
end