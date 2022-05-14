function L = mylaplacian1D(N,p,BC)
% L = mylaplacian1D(N,p,BC) returns the discretized sparse 1D laplacian 
% of order 2*p+2 and boundary conditions BC, that can be:
% 'PBC'   Periodic BC, the default;
% 'DBC0'  all zeros Dirichlet BC) namely f(j)=0 for j<=0 and j>=N+1;
% 'DBC'   Dirichlet BC: f(-j)=-f(j) and f(N+1+j)=-f(N+1-j) for j>=0;
% 'DSTI'  DBC through discrete sine transform of type I;
% 'NBC'   Neumann BC: f(1-j)=f(j) and f(N+j)=f(N-j+1) for j>0;
% 'DCTII' NBC through discrete cosine transform of type II,
% If p = inf, the dense Laplacian with exact k^2 spectrum is returned. In
% this case 'DBC0' is treated as 'DBC'. 

narginchk(2,3)
if nargin < 3
    BC = 'PBC'; 
end

if isfinite(p)
    if N < 2*p+1
        error('required matrix too small for the order');
    end
    [~,cp,kh2] = mymultidiagL(p);
    switch BC 
        case 'PBC'
            cp=2*double(cp);
            L=cp(1)*diag(ones(N,1));
            for i=2:length(cp)
                L=L+cp(i)*diag(ones(N-i+1,1),i-1)+cp(i)*diag(ones(i-1,1),N+1-i);
            end
        case 'DBC0' 
            L = zeros(N,N+p);
            cp = 2*double(cp);
            for j =1:N
                L(j,j:j+p) = cp;
            end
            L = L(:,1:N);
        case 'DBC'
            L = zeros(N,N+p);
            cp = 2*double(cp);
            for j =1:N
                L(j,j:j+p) = cp;
            end
            for j = 1:floor(p/2)
                for k = p-j:-1:j+1
                    L(j,k) = L(j,k)-sum(L(j,k+2*j:2*j:p+j+1));
                end
                L(j,j) = L(j,j)-sum(L(j,3*j:2*j:p+j))/2;
                L(N+1-j:-1:N-p+j+1,N+1-j) = L(j,j:p-j);
            end               
            L = L(:,1:N);
        case 'DSTI'
            k = (pi/(N+1))*(1:N)';
            L = idst(kh2(k).*dst(eye(N)));
        case 'NBC'
            L = zeros(N,N+p);
            cp = 2*double(cp);
            for j =1:N
                L(j,j:j+p) = cp;
            end
            L(1,1) = L(1,1)+cp(2)/2;
            L(1,2:p) = L(1,2:p)+cp(3:end);
            L(N:-1:N-p+1,N) = L(1,1:p);
            for j = 2:ceil(p/2)
                for k = j+1:p-j+1
                    L(j,k) = (-1).^(0:2*(j-1))*L(1,k-j+1:k+j-1)';
                end
                L(j,j) = L(1,1)+(-1).^(1:2*(j-1))*L(1,2:2*(j-1)+1)'/2;
                L(N+1-j:-1:N-p+j,N+1-j) = L(j,j:p+1-j);
            end
            L = L(:,1:N);
        case 'DCTII'
            k = (pi/N)*(0:N-1)';
            L = idct(double(kh2(k)).*dct(eye(N)));
        otherwise
            error('unknown type of boundary conditions');
    end
    L = real(L+L')/2;
    L(abs(L)<1e-10) = 0;  
    L = sparse(L);
else
   switch BC    
        case 'PBC'
            n = floor(N/2);
            nn = floor((N-1)/2);
            k = (2*pi/N)*(-n:nn)';
            L = isfft(diag(k.^2)*sfft(eye(N)));        
        case 'DBC0'
            error('order must be finite, use BC = ''DSTI'' instead of ''DBC0''');
        case 'DBC'
            error('order must be finite, use BC = ''DSTI'' instead of ''DBC''');
        case 'DSTI'
            k = (pi/(N+1))*(1:N)';
            L = idst(diag(k.^2)*dst(eye(N)));
        case 'NBC'
            error('order must be finite, use BC = ''DCTII'' instead of ''NBC''');
        case 'DCTII'
            k = (pi/N)*(0:N-1)';
            L = idct(diag(k.^2)*dct(eye(N)));
        otherwise
            error('unknown type of boundary conditions');
   end
   L = real(L+L')/2;
   L(abs(L)<1e-10) = 0;     
end


