fluxes = zeros(100, 6);
meanPerms = fluxes;
Ps = linspace(0, 10, 100);
Ms = [0 1 10 100 1000];

fluxes(:, 1) = Ps';
meanPerms(:, 1) = Ps';

for i = 1 : 5
    for j = 1 : 100
        [ph, a, q] = base_state(10, 2/3, Ms(i), Ps(j));
        fluxes(j, i+1) = -q;

        k = @(x)(ph(x).^(-2/3));
        meanPerms(j, i+1) = integral(k, 0, a)/a;
    end

    plot(Ps, fluxes(:, i+1));
end

writematrix(fluxes, 'output/fig4_q.csv');
writematrix(meanPerms, 'output/fig4_kBar.csv');