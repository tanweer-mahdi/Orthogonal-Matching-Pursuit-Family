function y = OMP_MMV(ZC,YY,sparsity)

%% Getting started
% ZC = Sensing matrix/Dictionary
% YY = Received data (Multiple Measurement Vector, also known as MMV)
% corrputed by AWGN noise
% sparsity = Prior sparsity level (number of non-zero element in sparse
% vector)
% Initialization
R = YY; %residual initialization
A = []; %active set initalization
[L B] = size(YY);
t = 1;
while(1)
    cor = ZC'*R;
    [v ii] = max(sum(abs(cor),2));
    A = [A ii]; % updating active indices
    phi = ZC(:,A);
    mod = phi'*phi + 1e-8*eye(size(phi'*phi,2)); % Tikhonov regularization
    P = eye(L) - phi*inv(mod)*phi'; % Null space projection matrix
    R = P*YY; % Updating residual by projecting to Null space 
    if length(A) == sparsity
        break;
    end
    t = t+1;
end

y = A;
end