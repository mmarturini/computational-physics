function draw_void(s)
%draw_void viene chiamata nel live script rw
%riceve come input una struttura dati del tipo output del programma
%randomwalk2
%viene utilizzata per plottare in funzione del tempo la differenza tra
%l'area del quadrato di lato 2sqrt(t) e l'area media del cerchio che ha
%come raggio la distanza dall'origine del random walker


N = s.in.number_of_particles;
n_step = size(s.positions,2)/2;
s.positions = reshape(s.positions,2*N,n_step);
dt = s.in.time_step;
t = 0:dt:(n_step-1)*dt;
A_circle = pi*(s.positions(1:N,:).^2 + s.positions(N+1:end,:).^2);  
A_mean = mean(A_circle,1);
A_void = 4*t - A_mean;
plot(t,A_void)
p = polyfit(t,A_void,1);
fprintf('Teoric values:\nm = %1.5f\tq = 0\n',4-2*pi*(s.in.sigma^2)); %confronto dei valori teorici con i valori ottenuti dall'interpolazione
fprintf('Interpolated values:\nm = %1.5f\tq = %1.5f\n',p(1),p(2));


