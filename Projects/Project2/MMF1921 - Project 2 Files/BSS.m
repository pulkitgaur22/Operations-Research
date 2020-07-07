function  [mu, Q, R_squre_adj] = BSS(returns, factRet, lambda, K)
    
    % Use this function for the BSS model. Note that you will not use 
    % lambda in this model (lambda is for LASSO).
    %
    % You should use an optimizer to solve this problem. Be sure to comment 
    % on your code to (briefly) explain your procedure.
    
 
    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    
    % r=X*B
    X=[ones(size(factRet,1),1) factRet];
    n=size(X,2); % Number of factors, including constant
    
    % min(||r-X*B||^2) -> min(-2*r'*X*B + B'*(X'*X)*B), Q = X'*X
    Q=X.'*X;
   
    % Introduce yi to constrain the number of non-zero bi
    Q=[Q zeros(n); zeros(n,2*n)];
    
    % Empty matrix to receive model-based stock expected return
    returnsModel=returns*0;
    % Empty matrix to receive betas
    betas=zeros(size(returns,2),n);
    % Empty matrix to receive residual variance
    D = zeros(size(returns,2),1);
    
    for i=1:size(returns,2) % Number of stocks
        rTemp=returns(:,i); % Stock under operation

        varTypes=[repmat('C',n,1);repmat('B',n,1)]; % n bi, n yi
        
        % -100*yi <= bi <= 100*yi
        A=[-eye(n),-100*eye(n);eye(n),-100*eye(n)];
        b=zeros(2*n,1);
        
        % Sum of yi equals to K
        Aeq=[zeros(1,n),ones(1,n)];
        beq=K;
        
        % Unbounded bi, [0,1] for yi
        lb=[-ones(1,n)*100,zeros(1,n)];
        ub=[ones(1,n)*100,ones(1,n)];

        % Append '_c' to cont var. names
        %namesCont = cellfun(@(c)[c '_c'], tickers(1:n), 'uni', false); 
        % Append '_b' to binary var. names
        %namesBin = cellfun(@(c)[c '_b'], tickers(1:n), 'uni', false); 

        % Combine both name vectors
        %names = [namesCont namesBin];

        clear model;
        %model.varnames = names; % Assign the variable names
        
        % Gurobi accepts an objective function of the following form:
        % f(x) = (1/2) x' Q x + obj' x 
        % Recall our objective: min(-2*r'*X*B + B'*Q*B)
        model.Q = sparse(Q);
        model.obj = [-2*rTemp.'*X,zeros(1,n)];
        
        % Gurobi only accepts a single A matrix
        % And use sense attribute to identify ineq/eq constrains
        model.A = [sparse(A); sparse(Aeq)];
        model.rhs = [b; beq];
        model.sense = [ repmat('<', (2*n), 1) ; '=' ];
        
        % Define the variable type (continuous and binary)
        model.vtype = varTypes;

        
        clear params;
        params.TimeLimit = 100;
        params.OutputFlag = 0; % Otherwise it output each iteration

        results = gurobi(model,params); % Optimize
        returnsModel(:,i)=X*results.x(1:9); % Fill in the expeted return
        betas(i,:)=results.x(1:9);
        D(i)=var(X*results.x(1:9)-rTemp);
    end
    
    V = betas(2:9,:);
    y = returns;
    %D=diag((1/(size(y,1)-1-size(V,1)))*D);
    D=diag(D)*(size(returns,1)-1)/(size(returns,1)-K);
    % Sum of squres

    SS_tot = sum((y - mean(y,1)).^2,1);
    SS_reg = sum((X*betas' - mean(y,1)).^2,1);
    SS_res = sum((y - X*betas').^2,1);
    R_squre = 1 - SS_res./SS_tot;
    
    R_squre_adj = 1 - (1-R_squre)*(size(y,1)-1)/(size(y,1)-1-size(V,1));
    mu = mean(returnsModel,1)'; % n x 1 vector of asset exp. returns
    Q = betas*cov(X)*betas.'+D; % n x n asset covariance matrix
    %Q = cov(returnsModel); % n x n asset covariance matrix

    %----------------------------------------------------------------------
    
end
