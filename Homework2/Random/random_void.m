function random_void(s)

%random_void riceve come input una struttura del tipo output del programma
%randomwalk2.m, contenente sia le informazioni relative alla creazione dei
%random walker sia le posizioni di ogni singolo random walker per ogni
%istante.

%Nota: la seguenete funzione viene utilizzata 'a pezzi' nel live script
%Rw.mlx (che da solo svolge tutto il procedimento), la riportiamo qua per completezza.

%la struttura letta da random_void salva le posizioni incolonnate una di
%fianco all'altra: le prime due colonne sono rispettivamente la x e la y al primo istante,
%le seconde due al secondo istante ecc...

dim = length(s.in.plot.viewbox)/2;
if dim ~= 2
    error('Random walk must be in 2 dimensions');
end

% parameters
N = s.in.number_of_particles;
n_step = size(s.positions,2)/2;
s.positions = reshape(s.positions,2*N,n_step); %reshape delle posizioni, posiziona in una colonna x e y, con le y sotto le x
dt = s.in.time_step;
t = 0:dt:(n_step-1)*dt;
A_circle = pi*(s.positions(1:N,:).^2 + s.positions(N+1:end,:).^2); %calcola l'area raggiunta da ogni random walker ad ogni istante 
%come superficie del cerchio cui raggio è distanza da origine
A_mean = mean(A_circle,1); %calcola media area raggiunta ad ogni istante
A_void = 4*t - A_mean;
diff2 = (A_circle - A_mean).^2;
sigma = sqrt(sum(diff2,1)/(N+1)); %calcola la sigma delle aree ad ogni istante
bins = linspace(0,5,30);

f1 = figure(1);
set(f1,'position',[100 120 700 600]);
clf
% plotta istoramma: mostra ad ogni istante la distribuzione delle aree 
a1 = subplot(2,2,1);
hi = histogram(A_circle(:,1),bins,'Normalization','pdf');
title('Area')
xlabel('Area covered');
ylabel('pdf');
ttl1 = title(sprintf('t = %-5.2f',0));
a1.YLim = [0 2]; 
pause(1)
for n = 2:n_step 
    hi.Data = A_circle(:,n); %aggiorna istogramma     
    if mod((n-1)*dt,0.01)
        ttl1.String = sprintf('t = %-5.2f',(n-1)*dt);
    end
    drawnow
    pause(0.001)
end

% plotta la sigma delle aree in funzione del tempo
subplot(2,2,2)
plot(t,sigma)
title('Standard deviation')
xlabel('Time');
ylabel('Sigma');
p = polyfit(t,sigma,1);
fprintf('Teoric values:\nm = %1.5f\tq = 0\n',2*pi*s.in.sigma^2);
fprintf('Interpolated values:\nm = %1.5f\tq = %1.5f\n',p(1),p(2));

% plotta ad ogni istante differenza tra area quadrato di raggio 2sqrt(t) e area media
% raggiunta da random walker
subplot(2,2,3)
plot(t,A_void)
title('Void area')
xlabel('Time');
ylabel('Difference');
p = polyfit(t,A_void,1);
fprintf('Teoric values:\nD = %1.5f\tq = 0\n',(s.in.sigma^2)/2);
fprintf('Interpolated values:\nD = %1.5f\tq = %1.5f\n',(4-p(1))/(4*pi),p(2));

% plotta il percorso dei primi tre random walkers
subplot(2,2,4)
plot(s.positions(1,:),s.positions(N+1,:))
hold on
plot(s.positions(2,:),s.positions(N+2,:))
plot(s.positions(3,:),s.positions(N+3,:))
title('')
xlabel('x');
ylabel('y');

% istogramma per singole storie (prime 20 particelle),a differenza del
% primo istogramma, ogni 'fotogramma' rappresenta la storia totale di un
% singolo random walker.

f2 = figure(2);
set(f2,'position',[1000 300 420 420]);
clf
a2 = axes;
hi = histogram(A_circle(1,:),bins,'Normalization','pdf');
title('Area singolo rw')
xlabel('Area covered');
ylabel('pdf');
a2.YLim = [0 2]; 
for n = 2:20
    i = randi(N);
    hi.Data = A_circle(i,:); %aggiorna istogramma
    drawnow
    pause(0.5)
end

end

