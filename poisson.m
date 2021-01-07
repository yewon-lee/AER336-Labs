% Copyright (c) 2016 Masayuki Yano, University of Toronto
classdef poisson
    methods (Access = public, Static = true)
        function [A,b] = getmatvec(m,ccflag)
            if (nargin < 2) ccflag = true; end
            [xq,wq] = poisson.gauss(2*m);
            [phi,phix] = poisson.lagrange(m+1,xq,ccflag);
            phi = phi(:,2:end-1);
            phix = phix(:,2:end-1);
            
            A1 = phix'*bsxfun(@times,wq,phix) + phi'*bsxfun(@times,wq,phi);
            b1 = phi'*(wq.*(xq > 0.5));
            A = kron(A1,A1);
            b = kron(b1,b1);
        end
        function vizsoln(u,ccflag)
            if (nargin < 2) ccflag = true; end
            m = sqrt(length(u));
            uu = zeros(m+2,m+2);
            uu(2:end-1,2:end-1) = reshape(u,[m,m]);
            xplot = linspace(0,1,30)';
            phi = poisson.lagrange(m+1,xplot,ccflag);
            uuplot = phi*uu*phi';
            surf(xplot,xplot,uuplot);
            %zlim([0,0.016]);
        end
    end
    methods (Access = private, Static = true)
        function [xq,wq] = gauss(m)
            beta = 0.5./sqrt(1-(2*(1:m-1)).^(-2));
            T = diag(beta,1) + diag(beta,-1);
            [~,D] = eig(T);
            xq = sort(diag(D));
            xq = 0.5*(xq+1);
            P0 = poisson.legendre(m-1,xq);
            A = inv(P0);
            wq = A(1,:).';
        end        
        function [phi,phix] = lagrange(p,x,ccflag)            
            xp = linspace(0,1,p+1)';
            if (ccflag == true)
                xp = 0.5*(1-cos(pi*xp));
            end
            psi0 = poisson.legendre(p,xp);
            [psi,psix] = poisson.legendre(p,x);
            phi = psi/psi0;
            phix = psix/psi0;
        end
        function [P,dP] = legendre(p,x)
            x = 2.0*x - 1.0; % x is over [0,1]; rescale to [-1,1]
            np = size(x,1);
            P = zeros(np,p+1);
            P(:,1) = 1.0;
            if (p > 0)
                P(:,2) = x;
            end
            for n = 1:p-1
                P(:,n+2) = (2*n+1)./(n+1)*x.*P(:,n+1) - n./(n+1)*P(:,n);
            end
            if (nargout > 1) % compute gradient
                ind = (abs(x) ~= 1);
                if any(ind)
                    dP(ind,2:p+1) = bsxfun(@rdivide,bsxfun(@times,1:p,(bsxfun(@times,x(ind),P(ind,2:p+1)) - P(ind,1:p))),x(ind).*x(ind)-1);
                end
                ind = (abs(x) == 1);
                if any(ind)
                    dP(ind,2:p+1) = repmat(((-1).^(0:p-1)).*(1:p).*(2:p+1)./2,[sum(ind),1]);
                end
                dP = 2.0*dP; % account for the [0,1] to [-1,1] scaling
            end
        end
    end
end
