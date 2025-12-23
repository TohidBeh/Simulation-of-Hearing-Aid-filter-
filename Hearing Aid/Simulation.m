clear all;
clc

%z = -10:0.01:10-0.01;
alpha = -0.657; 
beta = 1;
mu = -(1 + alpha)/(1 + beta);
syms z
r_z = mu * (1 + beta * z.^(-1)) ./ (1 + alpha * z.^(-1));
[Z , p] = tf2zp([mu mu*beta] , [1 alpha]);
figure(1)
zplane(Z , p);

r_n = iztrans(r_z);
figure(2)
% r_n_new = subs(r_n , n , 0:10);
hold on
fplot(r_n , 'r')
title("r(t) : ");
% t = 0:0.01:10;
% r_t = spline(-10:10 , t , r_n_new);
% plot(r_n_new)



%% part2:
K =  ((1 + mu) + (alpha + mu*beta) * z.^(-1)) ./ ((1 - mu) + (alpha + beta/mu) * z.^(-1));
[Z , p] = tf2zp([(1 + mu) (alpha + mu*beta)] , [(1 - mu) (alpha + beta/mu)]);
figure(3)
zplane(Z , p);
K_ap = ((alpha + beta/mu)/(1-mu) + z.^(-1)) ./ ((1 + (alpha + beta/mu)/(1 - mu)) * z.^(-1));
K_min = (1 + mu)/(1 - mu) * (1 + (alpha + mu * beta)/(1 + mu) * z.^(-1)) ./ ((alpha + beta/mu)/(1 - mu) + z.^(-1));
syms w;
K_ap = subs(K_ap , z , exp(1i*w));
K_min = subs(K_min , z , exp(1i*w));
group_delay_ap = -diff(angle(K_ap)) ;
group_delay_min = -diff(angle(K_min));

figure(4)
subplot(3,1,1)
fplot(20*log10(abs(K_ap)))
title("Absolute of all_pass");

subplot(3,1,2)
fplot(w , angle(K_ap))
title("Phase of all_pass");

subplot(3,1,3)
fplot(group_delay_ap)
title("Group delay of all_pass");

figure(5)
subplot(3,1,1)
fplot(abs(K_min))
title("Absolute of minimum_phase");

subplot(3,1,2)
fplot(angle(K_min))
title("Phase of minimum_phase");

subplot(3,1,3)
fplot(group_delay_min)
title("Group delay of minimum_phase");

