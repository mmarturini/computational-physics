function Ydot = lorenz_plusJ(~,Y,params)
% Output data must be a 12-element column vector
%mi da come output il vettore che rappresenta l'eq differenziale di lorenz
%per le prime tre componenti, mentre le rimanenti 9 componenti
%rappresentano le entrate della matrice che risolve il problema
%linearizzato Ydot = J*Y (vedi linearizz in intorno punto eq)
%
%  Lorenz equation 
%
%               dx/dt= - beta*x + y*z
%               dy/dt = sigma*(z - y)
%               dz/dt = -x*y + rho*y - z 
%
% plus Jacobian linearization 
% Output data must be a 12-element column vector


% Values of parameters
BETA = params(1);
RHO = params(2);
SIGMA = params(3);
     
x = Y(1); y = Y(2); z = Y(3);
Ydot = zeros(12,1);

%Lorenz equation
Ydot(1) = - BETA*x + y*z;
Ydot(2) = SIGMA*(z - y);
Ydot(3) = -x*y + RHO*y - z;

Y = reshape(Y(4:12),3,3);

%Jacobian of linearized system
J = [-BETA z y; 0 -SIGMA SIGMA; -y RHO-x -1];     
  
%Variational equation   
Ydot(4:12) = J*Y;

%J dipende dal tempo, quindi è un'evoluzione che dipende dal tempo

%se non dipendesse dal tempo avrei le solite soluzioni che ho per un
%sistema linearizzato in un intorno di un punto fisso.
%nel nostro caso invece non sono costanti ma dipendono dal tempo, e mylyap
%permette di calcolare gli exp di ljapunov quando  ho dipendenza temporale