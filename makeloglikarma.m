function ll = makeloglikarma(y, y1, n, k, m)
ll = @loglik;

    function [yy, J] = loglik(ini)
        
        
        alpha = ini(1,1);
        phi = ini(1,2:4);
        theta = ini(1,5:7);
        
        ynew = log(y);
        
        mu = zeros(n,k);
        eta = mu;
        error = mu;
        
        
        for i = (m+1):n
            
            for j = (m+1):k
                
                eta(i,j)  = alpha + phi(1,1)*ynew(i,j-1) + phi(1,2)*ynew(i-1,j) + ...
                    phi(1,3)*ynew(i-1,j-1) + ...
                    theta(1,1)*error(i,j-1) + theta(1,2)*error(i-1,j) + ...
                    theta(1,3)*error(i-1,j-1);
                
                mu(i,j) = exp(eta(i,j));
                
                error(i,j) = ynew(i,j)-eta(i,j);
                
                
            end
            
        end
        
        
        mu = (reshape(mu((m+1):n,(m+1):k)',1,(n-1)*(k-1)))';
        
        yy = -sum(sum(log(pi/2)+log(y1)-log(2.*mu.^2)-(pi.*y1.^2)./(4.*(mu.^2))));
        
        dmu = ((pi.*(y1.^2))./(2.*(mu.^3))-(2)./(mu));
        
        mT = diag(mu);
        
        
        
        deta_dalpha = zeros(n,k);
        deta_dphi1 = zeros(n,k);
        deta_dphi2 = zeros(n,k);
        deta_dphi3 = zeros(n,k);
        deta_dtheta1 = zeros(n,k);
        deta_dtheta2 = zeros(n,k);
        deta_dtheta3 = zeros(n,k);
        
        for i = (m+1):n
            for j = (m+1):n
                
                deta_dalpha(i,j) = 1 - (theta(1,1)*deta_dalpha(i,j-1) + ...
                    theta(1,2)*deta_dalpha(i-1,j) + ...
                    theta(1,3)*deta_dalpha(i-1,j-1));
                
                deta_dphi1(i,j) = ynew(i,j-1) - (theta(1,1)*deta_dphi1(i,j-1) + ...
                    theta(1,2)*deta_dphi1(i-1,j) + ...
                    theta(1,3)*deta_dphi1(i-1,j-1));
                
                deta_dphi2(i,j) =  ynew(i-1,j) - (theta(1,1)*deta_dphi2(i,j-1) + ...
                    theta(1,2)*deta_dphi2(i-1,j) + ...
                    theta(1,3)*deta_dphi2(i-1,j-1));
                
                
                deta_dphi3(i,j) =  ynew(i-1,j-1) - (theta(1,1)*deta_dphi3(i,j-1) + ...
                    theta(1,2)*deta_dphi3(i-1,j) + ...
                    theta(1,3)*deta_dphi3(i-1,j-1));
                
                deta_dtheta1(i,j) =  error(i,j-1) - (theta(1,1)*deta_dtheta1(i,j-1) + ...
                    theta(1,2)*deta_dtheta1(i-1,j) + ...
                    theta(1,3)*deta_dtheta1(i-1,j-1));
                
                deta_dtheta2(i,j) = error(i-1,j) - (theta(1,1)*deta_dtheta2(i,j-1) + ...
                    theta(1,2)*deta_dtheta2(i-1,j) + ...
                    theta(1,3)*deta_dtheta2(i-1,j-1));
                
                
                deta_dtheta3(i,j) = error(i-1,j-1) - (theta(1,1)*deta_dtheta3(i,j-1) + ...
                    theta(1,2)*deta_dtheta3(i-1,j) + ...
                    theta(1,3)*deta_dtheta3(i-1,j-1));
                
            end
        end
        
        
        a = reshape(deta_dalpha((m+1):n,(m+1):k)',1,(n-1)*(k-1));
        
        rP1 = reshape(deta_dphi1((m+1):n,(m+1):k)', 1,(n-1)*(k-1));
        rP2 = reshape(deta_dphi2((m+1):n,(m+1):k)', 1,(n-1)*(k-1));
        rP3 = reshape(deta_dphi3((m+1):n,(m+1):k)',1,(n-1)*(k-1));
        
        rR1 = reshape(deta_dtheta1((m+1):n,(m+1):k)', 1,(n-1)*(k-1));
        rR2 = reshape(deta_dtheta2((m+1):n,(m+1):k)', 1,(n-1)*(k-1));
        rR3 = reshape(deta_dtheta3((m+1):n,(m+1):k)', 1,(n-1)*(k-1));
        
        Ualpha = - a * mT * dmu;
        
        Uphi1 =   - rP1 * mT * dmu;
        Uphi2 =   - rP2 * mT * dmu;
        Uphi3 =   - rP3 * mT * dmu;
        
        Utheta1 =   - rR1 * mT * dmu;
        Utheta2 =   - rR2 * mT * dmu;
        Utheta3 =   - rR3 * mT * dmu;
        
        
        J = [Ualpha Uphi1 Uphi2 Uphi3 Utheta1 Utheta2 Utheta3];
        
        
    end


end