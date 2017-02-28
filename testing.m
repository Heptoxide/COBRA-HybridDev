%% Anthony Torres
% Testing

t = linspace(0,10);
mox = 0.5*exp(-(0.5.*t).^(1/3) - .3);

figure;
plot(t, mox)