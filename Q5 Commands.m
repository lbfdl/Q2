[un07, kn07, en07, tn07, pn07] = MDLJ(6, 0.8442, 0.7, 5150, 0.0222);
[un10, kn10, en10, tn10, pn10] = MDLJ(6, 0.8442, 1.0, 5150, 0.0222);
[un15, kn15, en15, tn15, pn15] = MDLJ(6, 0.8442, 1.5, 5150, 0.0222);
E07 = en07(150:end);
E10= en10(150:end);
E15= en15(150:end);
meanE07 = mean(E07);
meanE2_07 = mean(E07.^2);
Cv07_red = (meanE2_07 - meanE07^2) / 0.7;
% T* = 1.0
meanE10 = mean(E10);
meanE2_10 = mean(E10.^2);
Cv10_red = (meanE2_10 - meanE10^2) / 1.0;
% T* = 1.5
meanE15 = mean(E15);
meanE2_15 = mean(E15.^2);
Cv15_red = (meanE2_15 - meanE15^2) / 1.5;
epsilon_eV = 0.0031;
epsilon_J = epsilon_eV * 1.602e-19;  % eV to J
kB = 1.380649e-23; % J/K
Av_M = 3.353e-26; % (mole/kg) for SI
% Reduced sıcaklıklar
T_star = [0.7, 1.0, 1.5];
% Cv (reduced) değerlerin
Cv_reduced = [Cv07_red, Cv10_red, Cv15_red];
figure;
plot(T_star, Cv_reduced, '-o', 'LineWidth', 2);
xlabel('T*');
ylabel('Reduced Heat Capacity Cv* ');
title('Cv vs T in reduced units');
grid on;
% Real sıcaklıkları hesapla (K)
T_real = T_star .* (epsilon_J / kB);
% Real Cv (J/K.atom)
Cv_real = Cv_reduced.* kB/ Av_M;
% Plot: Cv (real) vs T (real)
figure;
plot(T_real, Cv_real, '-o', 'LineWidth', 2);
xlabel('Temperature (K)');
ylabel('Specific Heat Capacity Cv (J/K·kg)');
title('Cv vs T in Real Units');
grid on;

