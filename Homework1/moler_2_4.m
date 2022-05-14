function moler_2_4
% moler_2_4 utilizza 'bslashtx', una funzione che implementa
% un algorimtmo di risoluzione di sistemi lineari.
% in questo caso si risolve un circuito elettrico utilizzando due metodi
% differenti: trovando prima le correnti note le resistenze e la differenza
% di potenziale del source, e successivamente trovando le differenze di
% potenziale ai nodi conoscendo le conduttanze.
% Ricavando infine i voltaggi dalle correnti ottenute e viceversa
% confrontiamo i due risultati verificando l'equivalenza dei due metodi.

%si definisce il vettore delle correnti (valori scelti casualmente, range da 1 a 10)

r=randi(10,8,1);
g=1./r;
vs=5;
b=[0; 0; 0; vs];
c=[0;0;g(7)*vs;0];

%si definisce la matrice delle resistenze

R =[r(5)+r(1)+r(3)+r(8) -r(1) -r(3) -r(8);
    -r(1) r(4)+r(1)+r(2) -r(2) 0;
    -r(3) -r(2) r(2)+r(3)+r(6) -r(6);
    -r(8) 0 -r(6) r(7)+r(8)+r(6)];

%si ottiene il vettore delle correnti risolvendo il sistema lineare

i=bslashtx(R,b);

%si definisce la matrice delle conduttanze

G=[g(1)+g(2)+g(3) -g(1) -g(2) -g(3);
    -g(1) g(1)+g(4)+g(5) -g(4) 0;
    -g(2) -g(4) g(2)+g(4)+g(6)+g(7) -g(6);
    -g(3) 0 -g(6) g(3)+g(6)+g(8)];

%si ottiene il vettore dei voltaggi risolvendo il sistema lineare

v=bslashtx(G,c);

%si ricavano i voltaggi dalle correnti nel vettore 'i' e le correnti dai
%voltaggi nel vettore 'v'

V=[i(1)*(r(5)+r(1))-r(1)*i(2);
    r(5)*(i(1));
    vs-r(7)*i(4);
    r(8)*(-i(1)+i(4))];

I=[v(2)*g(5);
    (v(3)-v(2))*g(4);
    g(6)*(-v(3)+v(4))+v(4)*g(8)+v(2)*g(5);
    v(4)*g(8)+v(2)*g(5)];

% si osserva la differenza tra i voltaggi e le correnti 
% ottenuti con metodi differenti. La massima differenza ottenuta
% è dell'ordine di eps
formatSpec = 'La massima differenza per i voltaggi è %d\nLa massima differenza per le correnti è %d\n';
fprintf(formatSpec,max(abs(V-v)),max(abs(I-i)))






