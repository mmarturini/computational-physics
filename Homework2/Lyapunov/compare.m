function compare(par1,par2)
%Questa funzione è stata esplicitamente costruita per il secondo blocco del
%live script, e non ha carattere generale.
%Obiettivo della funzione è fornire una descrizione qualitativa della
%sensibilità alle condizioni iniziali del doppio pendolo.
%Fornisce un plot con l'evoluzione temporale di theta1 e un plot con quella
%di theta2; per ogni plot sono rappresentati due grafici corrispondenti
%ai due diversi angoli iniziali.
%caso 1
%struttura par1 del tipo generato da findlyap
%caso 2
%par1 vettore di angoli iniziali(esplicitamente 2);
%par2 vettore parametri masse e lunghezze;

switch nargin
    case 1
        swlyap = par1;
    case 2
        if length(par1)>2 
            error('Error: Too many input arguments')
        end
        swlyap = find_lyap(par1,par2);
end

clf
for j = 1:length(swlyap.theta)
    ax(1) = subplot(1,2,1);
    plot(swlyap.T,swlyap.out(j,1).path1)
    hold on
    axis square
    xlabel('t');
    ylabel('$\theta_1$','Interpreter','latex');
    
    ax(2) = subplot(1,2,2);
    plot(swlyap.T,swlyap.out(j,1).path2)
    hold on
    axis square
    xlabel('t');
    ylabel('$\theta_2$','Interpreter','latex'); 
    
    t(1) = title(ax(1),'Traiettoria $\theta_1(t)$');
    t(2) = title(ax(2),'Traiettoria $\theta_2(t)$');
    set(t,'Interpreter','latex');
end

%Stima tempo di Lyapunov
m = sort(swlyap.out(1).mean);%vettore con le medie dei coeff di Lyap del primo angolo
t = 1/m(end);                %valor medio maggiore dei coeff di Lyapunov
fprintf('Tempo di Lyapunov = %2.3f s',t);







