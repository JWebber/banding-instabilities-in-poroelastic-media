fluxes = zeros(100, 6);
meanPerms = fluxes;
Ps = logspace(0, 2, 100);
Ms = [0 1 10 100 1000];

fluxes(:, 1) = Ps';
meanPerms(:, 1) = Ps';

reference_fluxes = zeros(100, 11);
reference_fluxes(:, 1) = Ps';

j = 2;

for M = linspace(0.5, 2, 7)
    semilogx(Ps, M*Ps, 'LineWidth', 0.5, 'Color', 'k');
    hold on;
    reference_fluxes(:, j) = M*Ps;
    j = j+1;
end

writematrix(reference_fluxes, 'output/fig4_q_ref.csv');

for i = 1 : 5
    for j = 1 : 100
        [ph, a, q] = base_state(10, 2/3, Ms(i), Ps(j));
        fluxes(j, i+1) = -q;

        k = @(x)(ph(x).^(-2/3));
        meanPerms(j, i+1) = integral(k, 0, a)/a;
    end

    semilogx(Ps, fluxes(:, i+1));
    hold on;
end

writematrix(fluxes, 'output/fig4_q.csv');
writematrix(meanPerms, 'output/fig4_kBar.csv');