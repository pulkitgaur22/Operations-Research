function  [mu, Q, R_squre_adj] = OLS(returns, factRet, lambda, K)
    
    % Use this function to perform an OLS regression. Note that you will 
    % not use lambda or K in thismodel (lambda is for LASSO, and K is for
    % BSS).
 
    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    y = returns;
    
    % Add intercept
    x = [ones(size(factRet,1),1) factRet];
    
    % Algebra computation of alpha and betas
    Coefficients = ((x' * x) \ x') * y;
    
    % Extract betas only
    V = Coefficients(2:9,:);
    epsilon = y - x*Coefficients;
    %D = (1/(size(y,1)-1-size(V,1)))*diag(cov(epsilon));
    D = diag(cov(epsilon))*(size(y,1)-1)/(size(y,1)-size(V,1)-1);
    F = cov(factRet);
    
    % Sum of squres
    SS_tot = sum((y - mean(y,1)).^2,1);
    SS_reg = sum((x*Coefficients - mean(y,1)).^2,1);
    SS_res = sum((y - x*Coefficients).^2,1);
    R_squre = 1 - SS_res./SS_tot;
    
    mu =   Coefficients'*mean(x,1)';       % n x 1 vector of asset exp. returns
    Q  =   V'*F*V + diag(D);       % n x n asset covariance matrix

    R_squre_adj = 1 - (1-R_squre)*(size(y,1)-1)/(size(y,1)-1-size(V,1));
    %----------------------------------------------------------------------
    
end