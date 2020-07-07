function  [mu, Q, R_squre_adj] = LASSO(returns, factRet, lambda, K)
    
    % Use this function for the LASSO model. Note that you will not use K 
    % in this model (K is for BSS).
    %
    % You should use an optimizer to solve this problem. Be sure to comment 
    % on your code to (briefly) explain your procedure.
    
    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    lambdaList = 0.001:0.001:1;
    estimates= [];
    %Through this loop, we will iterate over all the 20 stocks and find the
    %corresponding weights for each stock's LASSO model
    for i =1:size(returns,2)
        for jj=1:size(lambdaList,2)
            lambda=lambdaList(jj);
        %subset the returns for ith stock.
            y= returns(:,i);
            
            %pre-calculating the number of factors and number of observations
            %for the factors matrix
            totalFactors = size(factRet,2);
            totalObs= size(factRet,1);
            
            %adding a columns of 1's to the factRet to account for alpha
            X= [ones(totalObs,1) factRet];
            
            %As we added auxiliary variables, adding columns for each of the
            %auxiliary variable, each factor has one auxiliary variable
            updatedfactRet= [X zeros(totalObs,totalFactors+1)];
            
            %As per matlab's quadprog, we calculate Q
            %Q = 2*X'*X
            Q=  2*(updatedfactRet')* updatedfactRet;
            
            %These are the set of 18 constraints
            %For each factor there is a set of 2 constraints,
            % y <=  B
            % y <= -B
            
            A= [eye(totalFactors+1)  -eye(totalFactors+1);
                -eye(totalFactors+1)  -eye(totalFactors+1)];
            
            %RHS vector of the 18 constraints
            b= zeros(size(A,1),1);
            
            %Computing cost vector now for the linear term
            %[ Bi's(9*1)  Aux. Variables(9*1) ]
            
            %calculating the cost of auxiliary variables
            %[ 0 0 0 0 0 0 0 0 0 lambda lambda ... lambda]
            auxiliaryVariables= [zeros(totalFactors+1,1) ; lambda*ones(totalFactors+1,1)]';
            
            %The weights for the factor is -2*Ri*X', for auxiliary zero
            %Adding the two set of weights, to get the final cost 
            c= ((-2*y'*updatedfactRet) + auxiliaryVariables)' ;
            
            %Setting options
            options = optimoptions('quadprog','TolFun',1e-9,'Display','off');
            
            
            %Using MATLAB's Quadprog calculating the optimum coefficients
            x= quadprog(Q,c,A,b,[],[],[],[],[],options);
            x= round(x(1:9),6);
            
            if (nnz(x)<=5 && nnz(x)>=2)
                estimates= [estimates x];
                break
            end
           
        end
    end
    y = returns;
    V = estimates(2:9,:);
    epsilon = returns - (X*estimates);
    %D = (1/(size(y,1)-1-size(V,1)))*diag(cov(epsilon));
    D = diag(cov(epsilon))*(size(y,1)-1)/(size(y,1)-nnz(V(1,:))-1);
    F = cov(factRet);
    
    % Sum of squres
    SS_tot = sum((y - mean(y,1)).^2,1);
    SS_reg = sum((X*estimates - mean(y,1)).^2,1);
    SS_res = sum((y - X*estimates).^2,1);
    R_squre = 1 - SS_res./SS_tot;
    
    mu =   estimates'*mean(X,1)';       % n x 1 vector of asset exp. returns
    Q  =   V'*F*V + diag(D);       % n x n asset covariance matrix
    R_squre_adj = 1 - (1-R_squre)*(size(y,1)-1)/(size(y,1)-1-size(V,1));
    %----------------------------------------------------------------------
    
end
