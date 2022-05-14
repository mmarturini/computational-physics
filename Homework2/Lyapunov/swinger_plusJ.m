function Tdot = swinger_plusJ(~,T,par)

%il T di input rappresenta [th1 th2 th1punto th2punto] a cui poi aggiungo i
%risultati dell'azione del jacobiano sarebbe stato meglio Y come in
%lorenplusJ

% Equazione del doppio pendolo e linearizzazione tramite matrice jacobiana
% L'output è un vettore colonna di 20 elementi, come richiesto da
% swinger_lyapexp; i primi 4 sono eq diff del doppio pendolo, i 16 sono
% invece sol del sist lineare

% Valori dei parametri
m1 = par(1);
m2 = par(2);    % masse in kilogrammi
L1 = par(3);
L2 = par(4);    % lunghezze in metri
    
u = T(1:4);
Tdot = zeros(20,1);

%dalla risoluzione della lagrangian del doppio pendolo ottengo una coppia
%di eq diff non lin e del secondo ordine per le eq del moto
%per riscriverle in forma lineare introduco il vettore u =[th1 th2 th1punto
%th2punto] di modo che se chiamo f la force function e M la matrice di massa
%le equazioni del moto diventano Mu.= f  (u. = u punto, cioè la forma
%vettorizzata che rappresenta le eq del moto)
%quindi per ottenere u. basta fare il backslash

% Costruzione della matrice di massa
c = cos(u(1)-u(2));
M = [1 0 0 0; 0 1 0 0; 0 0 (m1+m2)*L1 m2*L2*c; 0 0 m2*L1*c m2*L2];

% Campo di forze
g = 9.81;
s = sin(u(1)-u(2));
f = [u(3); u(4); -(m1+m2)*g*sin(u(1))-m2*L2*s*u(4)^2; 
    -m2*g*sin(u(2))+m2*L1*s*u(3)^2];
F = M\f;

% Equazione del doppio pendolo
Tdot(1:4) = F;

% Costruzione del jacobiano: è calcolato numericamente, perché la funzione
% simbolica jacobian richiede un tempo eccessivo
%da f=[f1 f2 f3 f4] calcolo la matrice jacobiana rispetto u1 u2 u3 u4
K = m1 + m2*s^2;
A = 2*m2*c^2/K;
B = 2*s/K;
J31 = 1/(K*L1)*(-g*(m1+m2)*cos(u(1))-m2*c*B*f(3)+s*f(4)*(1+A));
J32 = 1/(K*L1)*(g*m2*c*cos(u(2))+m2*c*B*f(3)-s*f(4)*(1+A));
J33 = -B*m2*c*u(3);
J34 = -B*m2*L2*u(4)/L1;
J41 = 1/(K*L2)*(s*f(3)*(1-A)+(m1+m2)*c*(g*cos(u(1))+B*f(4)));
J42 = 1/(K*L2)*(-s*f(3)*(1+A)+(m1+m2)*(-g*cos(u(2))+B*c*f(4)));
J43 = B*(m1+m2)*L1*u(3)/L2;
J44 = B*m2*c*u(4);
J = [0 0 1 0; 0 0 0 1; J31 J32 J33 J34; J41 J42 J43 J44];

%i 16 elementi dal quinto al 20esimo di T vengono reshapati in mat 4x4
T = reshape(T(5:20),4,4);

% Equazione variazionale
Tdot(5:20) = J*T;
