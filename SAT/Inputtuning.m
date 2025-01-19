gamma = [0.01, 0.3, 0.7, 1, 1.5, 2, 3];
c = [0, 6.4, 12.8, 25.6, 51.2]/100;
h = figure; hold on;
for i = 1:numel(gamma)
    V1 = (1+c).^gamma(i);
    V2 = (1-c).^gamma(i);
    plot(c, V1, '.-', 'LineWidth',2);
    %plot(c, V2, '--.');
end
xlabel('c''');
ylabel('V');
lg = legend({'0.1', '0.3', '0.7', '1', '1.5', '2', '3'});
title(lg, '\gamma');
title('V = (1+c'')^\gamma');


h = figure; hold on;
for i = 1:numel(gamma)
    V1 = 1+c.^gamma(i);
    V2 = 1-c.^gamma(i);
    plot(c, V1, '.-', 'LineWidth',2);
    %plot(c, V2, '--.');
end
xlabel('c''');
ylabel('V');
lg = legend({'0.1', '0.3', '0.7', '1', '1.5', '2', '3'});
title(lg, '\gamma');
title('V = 1+c''^\gamma');