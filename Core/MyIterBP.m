function [Xc Ec]=MyIterBP(Sig,N)
p=200;
D = randn(p,N);
D = D ./ repmat( sqrt(sum(D.^2,1)), [p 1] );
Y= D*Sig';
z = D'*Y;
lambda = std( z(:) )/100;

% conditioning of the dictionary
[a,s,b] = svd(D); 
mu = 1/max(diag(s))^2;

%%% noiseless inversion
options.niter = 1000;
options.lambda_min = 0;
c = D'*Y; c = c(:);
options.lambda_max = mad(z)/0.6745;
% options.thresh_type = 'hard';
options.thresh_type = 'soft';
[Xc,Ec] = perform_iterative_thresholding(D,Y(:,1),options);
Xc=Xc';
