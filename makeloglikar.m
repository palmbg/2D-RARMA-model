function ll = makeloglikar(y, y1, n, k, m, reg1, reg2, reg3)
ll = @loglik;

    function [yy, J] = loglik(ini)
        
        
        alpha = ini(1,1);
        phi = ini(1,2:max(size(ini)));
        
        P1 = reg1;
        P2 = reg2;
        P3 = reg3;
        
        ynew = log(y);
        
        mu = zeros(n,k);
        eta = mu;
        
        
        for i = (m+1):n
            
            for j = (m+1):k
                
                eta(i,j)  = alpha + (phi(1,1)*ynew(i,j-1) + phi(1,2)*ynew(i-1,j) + ...
                    phi(1,3)*ynew(i-1,j-1));
                mu(i,j) = exp(eta(i,j));
                
            end
            
        end
        
        mu = (reshape(mu((m+1):n,(m+1):k)',1,(n-1)*(k-1)))';
        
        
        yy = -sum(sum(log(pi/2)+log(y1)-log(2.*mu.^2)-(pi.*y1.^2)./(4.*(mu.^2))));
        
        dmu = ((pi.*(y1.^2))./(2.*(mu.^3))-(2)./(mu));
        
        mT = diag(mu);
        
        Ualpha = - sum(mT * dmu);
        Uphi1 =   - P1 * mT * dmu;
        Uphi2 =   - P2 * mT * dmu;
        Uphi3 =   - P3 * mT * dmu;
        
        
        J = [Ualpha Uphi1 Uphi2 Uphi3];
        
    end


end