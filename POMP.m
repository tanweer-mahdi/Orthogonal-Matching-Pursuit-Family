function [sup, resnorm] = POMP(D,y,ki,L)

%% This function computes Piecewise Orthogonal Matching Pursuit
% **Inputs**
% D = Sensing matrix
% y = Measurements
% ki = Sparsity in each block
% L = length of each block
% **Outputs**
% sup = Piecewise support

%% Initialization
Q1=mean(sum(abs(D).^2,1));
Q = length(y);
sup = [];
abandon = [];
res = y;
%ki = 1:(N/L);
%ki = ones(1,N/L);
imax = 500;
minres = 1e-3;
set = 1:size(D,2);
modset = set;
B0 = [];
adandon = [];

%% Performing POMP
for i=1:imax
    %norm(res,2)
    % adandoning the redundant columns from the main index
    modset = set(find(ismember(set,abandon)<1)); % modset = Modified set, the set of available columns
    % Computing dot products
    corx = abs(D(:,modset)'*res)/(sqrt(res'*res/Q));
    [val ind] = max(corx);
    % updating support
    sup = [sup modset(ind)];
    B0 = [B0 D(:,modset(ind))];
    % identifying corresponding block
    if rem(modset(ind),L) == 0
        block = modset(ind)/L;
    else
        block = fix(modset(ind)/L)+1; %block id
    end
    temp = (block-1)*L+1:block*L; %block indices
    % checking if sparsity in that block is maximum
    if sum(ismember(temp,sup))== ki(block)
        abandon = [abandon temp(find(ismember(temp,sup)<1))]; % updating abandoned column indices
    end
    % residual update
    res = (eye(length(res))-B0*pinv(B0))*y;
    if length(sup) == sum(ki) || norm(res,2)<1e-20
        break;
    end
end
resnorm = norm(res,2);