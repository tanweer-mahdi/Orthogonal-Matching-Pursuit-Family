function params = bsbl(ypn,pilot,lds,M,N,dc,pa)
params = struct();
%% Constructing the sensing matrix
PP = kron(pilot,eye(M));
%indicator = zeros(N,1);
%indicator(uset) = 1;
%chvec = [];
nz_indices = [];
ch_indices = [];
for i=1:N
    nz_indices = [nz_indices; find(abs(lds(:,i))>0) + (i-1)*M];
    ch_indices = [ch_indices; (find(abs(lds(:,i))>0))'];
    %chvec =[chvec; indicator(i)*h(:,i).*lds(:,i)];
end

PP = PP(:,nz_indices); % Sensing matrix constructed


%% BOMP
tempy = ypn;
res = tempy(:); % initialized residual
l = 0;
block_ind = [];
auset = [];
for i=1:N
    block_ind = [block_ind;(i-1)*dc + 1: i*dc];
end
temp_block_ind = [];
while l<ceil(N*pa)
    corr = zeros(N,1);
    for i=1:N
        temp = PP(:,(i-1)*dc + 1: i*dc);
        corr(i) = norm(temp'*res);
    end
    [~,ind] = max(corr);
    auset = [auset ind];
    temp_block_ind = [temp_block_ind block_ind(ind,:)];
    hh = PP(:,temp_block_ind)\tempy(:);
    res = tempy(:) - PP(:,temp_block_ind)*hh;
    l = l + 1;    
end
auset = sort(auset);
hh = PP(:,sort(temp_block_ind))\tempy(:);
%% Estimated channel matrix 
est_he = zeros(M,length(auset));
for i=1:length(auset)
    est_he(ch_indices(auset(i),:),i) = hh((i-1)*dc+1:i*dc);
end

params.auset = auset;
params.est_channels = est_he;