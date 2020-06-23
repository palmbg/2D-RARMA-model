function [coef, mu_hat, resid1, resid2, pvalue] = arfit(phi1, phi2, y, n, k)

% AR model

p = max(max(size(phi1), size(phi2)));
m = p;

ar1 = 1:1:max(size(phi1));
ar2 = 1:1:max(size(phi2));

ynew = log(y);


X = zeros(1, 2 + ar1 + ar2);
XX = zeros( 1, 2 + ar1 + ar2);


for i = (m+1):n
    
    for j = (m+1):k
        
        XX(1,:) = [ynew(i,j) ynew(i,j-1) ynew(i-1,j) ynew(i-1, j-1)];
        X = [X; XX];
        
    end
    
    
end

nx = size(X);

P = X(2:nx(1,1), :);

y1 = ones((n-m)*(k-m),1);

x1 = [y1 P(:,2:nx(1,2))];

ini = regress(P(:,1), x1);


loglik = makeloglikar(y, exp(P(:,1)), n, k, m, P(:,2)', P(:,3)', P(:,4)');
options = optimoptions('fminunc', 'SpecifyObjectiveGradient', true, ...
    'Display', 'none',  'HessUpdate', 'bfgs', 'MaxIterations', 400, ...
    'OptimalityTolerance', 1e-6, 'ObjectiveLimit', -1e20);
[opt2,~,exitflag,~, ~, ~] = fminunc(loglik, ini', options);

if (exitflag ~= 0)
    
    coef = opt2;
    
    alpha = coef(1,1);
    phi = coef(1,2:4);
    
    etahat = zeros(n,k);
    muhati = etahat;
    
    for i = (m+1):n
        
        for j = (m+1):k
            
            etahat(i,j) = alpha + phi(1,1)*ynew(i,j-1) + phi(1,2)*ynew(i-1,j) + ...
                phi(1,3)*ynew(i-1,j-1);
            muhati(i,j) = exp(etahat(i,j));
            
        end
        
    end
    
    mu_hat = muhati((m+1):n,(m+1):k);
    
    y_1 = y((m+1):n,(m+1):k);
    
    var_y = (mu_hat.^2)*(4/pi -1)';
    
    resid1 = (y_1-mu_hat)./sqrt(var_y);
    
    cumulativef = pr(y_1,mu_hat);
    
    resid2 = norminv(cumulativef);
    
    mu_hat2 = (reshape(mu_hat',1,(n-1)*(k-1)))';
    
    W = diag(((4)./(mu_hat2.^2)).*(mu_hat2.^2));
    
    vi = ones(1,(n-m)*(k-m));
    
    
    Kaa =   vi * W * vi';
    Kap1 =  P(:,2)' * W * vi' ;
    Kap2 =  P(:,3)' * W * vi' ;
    Kap3 =  P(:,4)' * W * vi' ;
    
    Kp1a = Kap1;
    Kp1p1 = P(:,2)' * W * P(:,2) ;
    Kp1p2 = P(:,2)' * W * P(:,3) ;
    Kp1p3 = P(:,2)' * W * P(:,4) ;
    
    Kp2a = Kap2;
    Kp2p1 = Kp1p2;
    Kp2p2 = P(:,3)' * W * P(:,3) ;
    Kp2p3 = P(:,3)' * W * P(:,4) ;
    
    Kp3a = Kap3;
    Kp3p1 = Kp1p3;
    Kp3p2 = Kp2p3;
    Kp3p3 = P(:,4)' * W * P(:,4) ;
   
    
    K1 = [Kaa Kap1 Kap2 Kap3 ];
    K2 = [Kp1a Kp1p1 Kp1p2 Kp1p3];
    K3 = [Kp2a Kp2p1 Kp2p2 Kp2p3];
    K4 = [Kp3a Kp3p1 Kp3p2 Kp3p3];
    
    K = [K1; K2; K3; K4];
    
    vcov = inv(chol(K)) * (inv(chol(K)))';
    
    stderror = sqrt(diag(vcov));
    
    
    zstat = abs(coef(1,2:4)/stderror(2:4,1)');

    pvalue = 2*(1 - normcdf(zstat));
    
    
    
    
end

end



