clear;
clc;

mex('dpotrf_.c', '-lmwblas','-lmwlapack');
mex('dpotrs_.c', '-lmwblas','-lmwlapack');

load('test_data');
Xtrain = Data.Xtrain;
Xtest = Data.Xtest;
Ytrain = Data.Ytrain;
Ytest = Data.Ytest;

Iters = length(Xtrain);
Tasks = size(Ytrain{1,1},2);
Dimen = size(Xtrain{1,1},2);

ind = load('ind.txt');
indt = ind;
for i = 1:Tasks-1
        indt = cat(1,indt,ind+Dimen*i);
end
ind = [0,indt'];

% Group Property
l = length(ind)-1;
for i = 1:l
    a = ind(i);
    b = ind(i+1);
    nn = b-a;
    temp = [a+1,b,sqrt(nn)];
    ind_t(:,i) = temp';
end
ind = ind_t;

% set opts
opts.Print = 0;
opts.MaxIter = 1000;
opts.Tol = 1e-3;
opts.rho = 10;
opts.Dimen = Dimen;
opts.ind = ind;
opts.Tasks = Tasks;

lambda1 = 1;
lambda2 = 1;
corthreshold = 0;

for i = 1:Iters
    fprintf('Iters,%d\n',i);
    Xtr = Xtrain{1,i};
    Xte = Xtest{1,i};
    Ytr = Ytrain{1,i};
    Yte = Ytest{1,i};
    
    [E, C] = taskNetwork(Ytr,corthreshold);
    
    h = size(E,1);      
    C = -C + eye(Tasks);
    C = C + sum(abs(C(:)))*eye(Tasks);
    D = C/h; 
        
    [W] = FLADMM(Xtr, Ytr, lambda1, lambda2, D, opts);
end
    