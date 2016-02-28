clear;
rng('default');
path(genpath('~/GitHub/cs181-practicals/'), path)

if false
    load('circular_sub.mat');
    ind = load('circular_sub_ind.txt');
    d = dataset('file', 'train.csv', 'ReadObsNames', true, 'ReadVarNames', true, 'Delimiter', ',');

    x = m; clear m;
    y = d.gap(ind);
    
    % Remove features that are always zero or almost always 1
    torm = [];
    torm = [torm, find(sum(x) == 0)];
    torm = [torm, find(sum(x)>0.95*size(x,1))];


    % Remove redundant predictors
    R = corr(x);
    corr_thresh = 0.90;
    for i = 1 : size(R,2)
        for j = i + 1 : size(R,2) 
            if R(i,j) > corr_thresh
                torm = [torm,j];
            end
        end
    end


    x(:,unique(torm)) = [];
    
    save('circular_sub_processed.mat', 'x', 'y', 'torm');
else
    load('circular_sub_processed.mat');
end

if false
    % Pick basis set of molecules
    [xs,mux,sigx] = standardize(x);
    [c,s,~,~,pexp] = pca(xs);
%     [~,bimax] = max(s);
%     [~,bimin] = min(s);
%     ncomp = 500;
%     bi = unique([bimax(1:ncomp), bimin(1:ncomp)]);
    
    bi_total = randsample(size(x,1), 11000);
    bi = bi_total(1:10000);
    rind = bi_total(10001:11000);


    save('basis_idx.mat', 'bi', 'rind', 'c', 's', 'pexp', 'xs', 'mux', 'sigx');
    return
else
    load('basis_idx.mat');
end

xrot = xs*c; %clear s;

% Training - local regression model
% r = corr(x,y);
% mask = abs(r) > 0.2;

K = @gaussian_kernel_diagonal; %(1 - x*x0'/size(x,2))/lambda; %exp( -(1 - x*x0'/size(x,2))/lambda ); % RBF of a hamming distance
lambda = 0.01:0.02:0.3;
x_basis = xrot(bi,:);
y_basis = y(bi);

scale = sum(std(xrot));
params{1} = K;
params{2} = []; %*0.1;
params{3} = false;
params{4} = false;
params{5} = 'local_average';

cvacc = zeros(length(lambda), 2);
for i = 1 : length(lambda)
    fprintf('%0.0f/%0.0f CV iterations complete.', i, length(lambda));
    params{2} = scale*lambda(i);
    model = kernelfit_train(x_basis,y_basis, params);
    yhat = kernelfit_predict(xrot(rind,:),model);
    
    cvacc(i,1) = sqrt(var(y(rind)-yhat,1));
    cvacc(i,2) = corr(yhat, y(rind));
end



 

% d = dataset('file', 'train.csv', 'ReadObsNames', true, 'ReadVarNames', true, 'Delimiter', ',');
% x = double(d(:,1:256));
% y = d.gap;
% 
% % Remove features that are always zero or almost always 1
% fm0 = (sum(x) ~= 0);
% x0 = x(:,fm0);
% x0(:,sum(x0)>0.95*1000000) = [];
% 
%