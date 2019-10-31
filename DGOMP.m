function y = DGOMP(ZC,YY,var,PFA)

%% Getting started
% ZC = Sensing matrix/Dictionary
% YY = Received data (Multiple Measurement Vector, also known as MMV)
% corrputed by AWGN noise
% var = Noise variance, assumed known as a priori
% PFA = Probability of False Alarm. A good starting point is PFA = 0.1
% Initialization
R = YY; %residual initialization
A = []; %active set initalization
[L B] = size(YY);
t = 1;
while(1)
    cor = ZC'*R;
    [v ii] = max(sum(abs(cor),2));
    A = [A ii];
    phi = ZC(:,A);
    mod = phi'*phi + 1e-8*eye(size(phi'*phi,2));
    P = eye(L) - phi*inv(mod)*phi';
    R = P*YY;
    ll = [eye(L-t) zeros(L-t,t)]*P;
    Z = ll*R;
    % Computing Threshold
    COV = ll*ll' + 1e-8*eye(L-t);
    T = zeros(1,B);
    for i=1:B
        T(i) = Z(:,i)'*inv(COV)*Z(:,i)/(var);
    end
    
    threshold = chi2inv((1-PFA)^(1/B),L-t);
    if sum(abs(T)<threshold) == B
        break;
    end
    t = t+1;
end

y = A;
end