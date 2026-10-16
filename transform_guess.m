function solinit = transform_guess(sol)
    p = sol.y(1, :);
    pP = sol.y(2, :);

    e = sol.y(3, :);
    eP = sol.y(4, :);

    P = -eP;
    PP = -sol.y(5, :);

    guessVals = [p; pP; 0*p; 0*pP; e; eP; P; PP; sol.y(6, :); sol.y(7, :)];

    solinit = bvpinit(sol.x, @(x)interp1(sol.x, guessVals', x)', [sol.parameters(1)]);
end