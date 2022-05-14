function draw_lyapexp(par1,par2)
%La funzione fornisce i grafici dell'evoluzione temporale dei 4 coefficienti
%di Lyapunov, della traettoria di theta2 e dell'evoluzione temporale della 
%somma dei coefficienti.
%caso nargin=0
%default di find_lyap
%caso nargin=1
%par1 struttura del tipo generata da find_lyap
%caso nargin=2
%par1 vettore contenente diversi possibili angoli iniziali (in radianti)
%       par1=[0.8 1 1.7 ... ]
%par2 matrice nx4 contenente diversi possibili valori di masse e lunghezze
%       par2=[m11 m12 L11 L12; m21 m22 L21 L22; ... ]
%Per ogni angolo fornisce una figura con 6 subplot, per ogni sublot sono 
%rappresentati grafici per diverse masse e lunghezze

switch nargin
    case 0
        swlyap = find_lyap;
    case 1
        swlyap = par1;
    case 2
        swlyap = find_lyap(par1,par2);
end
%alla fine della fiera swlyap è una struttura del tipo generato da findlyap

ax = zeros(6,1);
t = zeros(6,1);

for j = 1:length(swlyap.theta) %numero di figure pari al numero di angoli passati
    figure(j)
    clf
  
 %per ogni figura
    %i 4 subplot dell'andamento temp dei coeff di lyap
    for i = 1:4 
        ax(i) = subplot(2,3,i);
        plot(swlyap.T,swlyap.out(j,1).L(:,i))
        hold on
        for k = 2:size(swlyap.par,1)
            plot(swlyap.T,swlyap.out(j,k).L(:,i))
        end  
        axis square
        xlabel('t');
        ylabel('Re $\lambda$','Interpreter','latex');
    end   
 
    %subplot traettoria del secondo angolo
    ax(5) = subplot(2,3,5);
    plot(swlyap.T,swlyap.out(j,1).path2)
    hold on
    for k = 2:size(swlyap.par,1)
         plot(swlyap.T,swlyap.out(j,k).path2)
    end    
    axis square
    xlabel('t');
    ylabel('$\theta_2$','Interpreter','latex');
    
    %subplot evo temp somma dei coeff di lyap
    ax(6) = subplot(2,3,6);
    plot(swlyap.T,swlyap.out(j,1).sum)
    lg=legend('par[1]','Location','southeast');
    hold on
    for k = 2:size(swlyap.par,1)
        plot(swlyap.T,swlyap.out(j,k).sum,'DisplayName',['par['...
             num2str(k) ']'])
    end  
    if size(swlyap.par,1) == 1
        set (lg,'visible','off');
    end 
    axis square
    xlabel('t');
    ylim([-1 0.05]);
    
    % titoli
    t(1) = title(ax(1),'$\lambda_1$');
    t(2) = title(ax(2),'$\lambda_2$');
    t(3) = title(ax(3),'$\lambda_3$');
    t(4) = title(ax(4),'$\lambda_4$');
    t(5) = title(ax(5),'Traiettoria $\theta_2(t)$');
    t(6) = title(ax(6),'Somma dei $\lambda$');
    set(t,'Interpreter','latex');  
    
    str = sprintf('%s',sym(swlyap.theta(j)));
    sgtitle(['Coefficienti di Lyapunov e traiettorie per \theta = ' str]);
end

fprintf('\t\t m1\t\t m2\t\t L1\t\t L2\n');
for k = 1:size(swlyap.par,1)
    fprintf('par[%d]\t\t',k);
    fprintf('%2.1f\t\t',swlyap.par(k,:));
    fprintf('\r');
end
     
%esempio di chiamate:
%draw_lyapexp

%s = find_lyap([pi/2 pi/4],[1 1 1 1])
%draw_lyapexp(s)

%draw_lyapexp([1.7 1.78],[1 1 1 1;2 2 2 2])





