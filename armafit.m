function [coef, mu_hat, resid1, resid2, pvalues] = armafit(y, n, k)

% ARMA model

m = 1;

ynew = log(y);

X = zeros(1, 4);
XX = zeros( 1, 4);


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

ini = [regress(P(:,1), x1); 0 ;0 ;0];


loglik = makeloglikarma(y, exp(P(:,1)), n, k, m);
options = optimoptions('fminunc', 'SpecifyObjectiveGradient', true, ...
    'Display', 'none',  'HessUpdate', 'bfgs', 'MaxIterations', 400, ...
    'OptimalityTolerance', 1e-6, 'ObjectiveLimit', -1e20);
[opt2,~,exitflag,~, ~, ~] = fminunc(loglik, ini', options);



if (exitflag ~= 0)
    
    coef = opt2;
    
    alpha = coef(1,1);
    phi = coef(1,2:4);
    theta = coef(1,5:7);
    
    etahat = zeros(n,k);
    muhati = etahat;
    errorhat = etahat;
    
    for i = (m+1):n
        
        for j = (m+1):k
            
            etahat(i,j)  = alpha + phi(1,1)*ynew(i,j-1) + phi(1,2)*ynew(i-1,j) + ...
                phi(1,3)*ynew(i-1,j-1) + ...
                theta(1,1)*errorhat(i,j-1) + theta(1,2)*errorhat(i-1,j) + ...
                theta(1,3)*errorhat(i-1,j-1);
            
            muhati(i,j) = exp(etahat(i,j));
            
            errorhat(i,j) = ynew(i,j) -log(muhati(i,j));
            
        end
        
    end
    
    mu_hat = muhati((m+1):n,(m+1):k);
    
    y_1 = y((m+1):n,(m+1):k);
    
    var_y = (mu_hat.^2)*(4/pi -1)';
    
    resid1 = (y_1-mu_hat)./sqrt(var_y);
    
    cumulativef = pr(y_1,mu_hat);
    
    resid2 = norminv(cumulativef);
    
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
            
            deta_dtheta1(i,j) =  errorhat(i,j-1) - (theta(1,1)*deta_dtheta1(i,j-1) + ...
                theta(1,2)*deta_dtheta1(i-1,j) + ...
                theta(1,3)*deta_dtheta1(i-1,j-1));
            
            deta_dtheta2(i,j) = errorhat(i-1,j) - (theta(1,1)*deta_dtheta2(i,j-1) + ...
                theta(1,2)*deta_dtheta2(i-1,j) + ...
                theta(1,3)*deta_dtheta2(i-1,j-1));
            
            
            deta_dtheta3(i,j) = errorhat(i-1,j-1) - (theta(1,1)*deta_dtheta3(i,j-1) + ...
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
    
    
    mu_hat2 = (reshape(mu_hat',1,(n-1)*(k-1)))';
    
    W = diag(((4)./(mu_hat2.^2)).*(mu_hat2.^2));
    
    
    Kaa =   a * W * a';
    Kap1 =  rP1 * W * a' ;
    Kap2 =  rP2 * W * a' ;
    Kap3 =  rP3 * W * a' ;
    Kat1 =  rR1 * W * a' ;
    Kat2 =  rR2 * W * a' ;
    Kat3 =  rR3 * W * a' ;
    
    Kp1a = Kap1;
    Kp1p1 = rP1 * W * rP1' ;
    Kp1p2 = rP1 * W * rP2' ;
    Kp1p3 = rP1 * W * rP3' ;
    Kp1t1 = rP1 * W * rR1' ;
    Kp1t2 = rP1 * W * rR2' ;
    Kp1t3 = rP1 * W * rR3' ;
    
    Kp2a = Kap2;
    Kp2p1 = Kp1p2;
    Kp2p2 = rP2 * W * rP2' ;
    Kp2p3 = rP2 * W * rP3' ;
    Kp2t1 = rP2 * W * rR1' ;
    Kp2t2 = rP2 * W * rR2' ;
    Kp2t3 = rP2 * W * rR3' ;
    
   
    Kp3a = Kap3;
    Kp3p1 = Kp1p3;
    Kp3p2 = Kp2p3;
    Kp3p3 = rP3 * W * rP3' ;
    Kp3t1 = rP3 * W * rR1' ;
    Kp3t2 = rP3 * W * rR2' ;
    Kp3t3 = rP3 * W * rR3' ;
    
    Kt1a = Kat1;
    Kt1p1 = Kp1t1 ;
    Kt1p2 = Kp2t1 ;
    Kt1p3 = Kp3t1 ;
    Kt1t1 = rR1 * W * rR1' ;
    Kt1t2 = rR1 * W * rR2' ;
    Kt1t3 = rR1 * W * rR3' ;
    
    Kt2a = Kat2;
    Kt2p1 = Kp1t2 ;
    Kt2p2 = Kp2t2 ;
    Kt2p3 = Kp3t2 ;
    Kt2t1 = Kt1t2 ;
    Kt2t2 = rR2 * W * rR2' ;
    Kt2t3 = rR2 * W * rR3' ;
    
    Kt3a = Kat3;
    Kt3p1 = Kp1t3 ;
    Kt3p2 = Kp2t3 ;
    Kt3p3 = Kp3t3 ;
    Kt3t1 = Kt1t3 ;
    Kt3t2 = Kt2t3 ;
    Kt3t3 = rR3 * W * rR3' ;
   
    
    K1 = [Kaa Kap1 Kap2 Kap3 Kat1 Kat2 Kat3];
    K2 = [Kp1a Kp1p1 Kp1p2 Kp1p3 Kp1t1 Kp1t2 Kp1t3];
    K3 = [Kp2a Kp2p1 Kp2p2 Kp2p3 Kp2t1 Kp2t2 Kp2t3];
    K4 = [Kp3a Kp3p1 Kp3p2 Kp3p3 Kp3t1 Kp3t2 Kp3t3];
    K5 = [Kt1a Kt1p1 Kt1p2 Kt1p3 Kt1t1 Kt1t2 Kt1t3];
    K6 = [Kt2a Kt2p1 Kt2p2 Kt2p3 Kt2t1 Kt2t2 Kt2t3];
    K7 = [Kt3a Kt3p1 Kt3p2 Kt3p3 Kt3t1 Kt3t2 Kt3t3];
    
    K = [K1; K2; K3; K4; K5; K6; K7];
    
    vcov = inv(chol(K)) * (inv(chol(K)))';
    
    stderror = sqrt(diag(vcov));
    
    
    zstat = abs(coef(1,2:7)/stderror(2:7,1)');

    pvalues = 2*(1 - normcdf(zstat));
    
end

end



