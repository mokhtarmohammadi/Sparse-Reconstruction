function x=MyBP(Sig,N)
p=2000;
D = randn(p,N);
D = D ./ repmat( sqrt(sum(D.^2,1)), [p 1] );

f = D*Sig';
options.tol = 1e-3;  % tolerance
options.nbr_max_atoms = 200;
options.use_slow_code = 0;
x1 = perform_mp(D,f,options);
x=x1';