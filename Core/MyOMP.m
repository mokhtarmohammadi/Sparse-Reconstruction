function x=MyOMP(Sig,iii,N)
m=N-(length(iii));
p=2000;
D = randn(p,N);
D = D ./ repmat( sqrt(sum(D.^2,1)), [p 1] );
f = D*Sig';
% options.tol = 1e-3;  % tolerance
% options.nbr_max_atoms = m;
% options.use_slow_code = 0;
%  options.sparse_coder = 'omp';
% x = perform_omp(D,f,options);
x=OMPa(D,f,m,N);

x=x';