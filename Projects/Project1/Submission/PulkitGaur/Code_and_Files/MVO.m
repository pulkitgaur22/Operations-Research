function  x = MVO(mu, Q, targetRet)
    
    % Use this function to construct your MVO portfolio subject to the
    % target return, with short sales disallowed. 
    %
    % You may use quadprog, Gurobi, or any other optimizer you are familiar
    % with. Just be sure to comment on your code to (briefly) explain your
    % procedure.

    % Find the total number of assets
    warning('off');
    n = size(Q,1); 

    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    
    H = Q; % Objective function cost
    % Constraints
    Aeq = ones(1,n); % sum of weights is equal to 1, portfolio return meet the required return
    beq = 1;
    A = -mu';
    b = -targetRet;
    lb = zeros(n,1); % short sales are disallowed
    
    %Setting options
    options = optimoptions('quadprog','TolFun',1e-30,'Display','off');

    x = quadprog(H,[],A,b,Aeq,beq,lb,[],[],options);        % Optimal asset weights
    %----------------------------------------------------------------------
    
end