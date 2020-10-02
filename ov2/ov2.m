G1 = tf([1 -2 ], [1 8 15]);

figure(1);
pzmap(G1);

figure(2);
nyquist(G1);