function sparsity_aware(Sig,P,S_L)
currentDirectory = pwd;
[upperPath, ~, ~] = fileparts(currentDirectory);
addpath([upperPath '\Core'])
addpath([upperPath '\TFSA7'])

%%%%%% Error calculation on nonsparse uniform STFT %%%%%%
%%%
%%%%%% The code is for the algorithm published in the paper:
%%%%%% "On the reconstruction of nonsparse time-frequency signals with
%%%%%% sparsity constraint from a reduced set of samples"
%%%%%% Isidora Stankovic, Cornel Ioana, Milos Dakovic
%%%%%% Signal Processing, vol. 142, pp. 480?484, January 2018.
%%%
%%%%%% When using the code, please cite this paper.
%%%%%% Code updated: December 2018.
%%%
%%%%%% Functions needed: SM_calc.m
addpath('D:\tfsa_5-5\windows\win64_bin');
WN=64; N=128;
n=0:N-1+WN*1;

Sig1 = exp(1i*(2*pi*(0.05*n.^3)/(N*N))+1i*(2*pi*(0.05*n))+1i*(2*pi*(0.05*n.^2/N))); %300t»???150t
Sig2= exp(1i*(2*pi*(0.05*n.^3)/(N*N))+1i*(2*pi*(0.2*n))+1i*(2*pi*(0.05*n.^2/N)));
%Sig2= exp(1i*(-2*pi*(0.1*n.^3)/(N*N))+1i*(2*pi*(0.45*n))-1i*(2*pi*(0.05*n.^2/N)));

%Sig1 = exp(1i*(2*pi*(0.1*n))); %300t»???150t
%Sig2= exp(1i*(2*pi*(0.3*n)));
% n=0:N-1;



x_o=Sig1+Sig2;

x_p=[x_o ];
%x_p=[zeros(1,WN/2) x_o zeros(1,WN/2)];
x_s=x_p;

%  S_L=56;

  S_L=32; % Number of missing samples
  p = randperm(N+WN);
  
   %    p=[60:70 91:100 110:120];
    
 x_s(p(1:S_L))=0;
% x_s=[ x_s ];
NA=find(x_s~=0);

STFT=[]; STFT_R=[];
it=0;
Sig_r=zeros(1,length(x_p));

wind_step=32;

for tt=0:WN/wind_step:length(x_p)-WN-1    % Signal length (1024+256) with a step of 32
    it=it+1;
    
    
    K=24;                % Sparsity
 
    x_w=x_p(tt+1:tt+WN).'.*hanning(WN,'periodic');
    
    
    X=(fft(x_w));                    % FT of the original (nonsparse) signal
           
    
    STFT=[STFT, X];              % STFT of the original (nonsparse) signal
    
    
    NA_v=NA(and(NA<tt+WN,NA>tt))-tt; % Random positions of the available samples
    
    % Generate the measurement matrix:
    IDFT_mtx=conj(dftmtx(WN))/WN; 

    A=IDFT_mtx(NA_v,:);             % Partial IDFT matrix

    y=A*X;                          % Available samples

    
    % Reconstruction with iterative OMP:
    KB=[];  y0=y; Kstep=1;
    for iter=1:round(K/Kstep)
        
        if iter>1; y=y0-x_R; end                 % Remove the reconstructed component from the array
        
        X0=(A*WN)'*y;                            % Initial estimate
        [Xv,Kv]=sort(abs(X0));                   % Sort the components
        KB=union(KB,Kv(end-Kstep+1:end));        % Take the largest component
        A_K=A(:,KB);                             % Measurement matrix for the largest component
        X_R=zeros(WN,1);
        X_R(KB)=pinv(A_K)*y0;                    % FT reconstruction of the components
        x_R=A*X_R;                               % Inverse of X_R
    end
    
    STFT_R=[STFT_R, X_R]; % Reconstructed STFT
    Sig_r(tt+1:tt+WN)=Sig_r(tt+1:tt+WN)+ifft(X_R).'/(wind_step/2);
    
    % Statistical error per one window:        
    E_stat(it)=10*log10(sum(abs(X_R(KB)-X(KB)).^2));
    
end

%  [extr_Sig2]=extridge_mult_tffiltering_multiple_tfd_viterbi(x_s, 2, 5,3,x_o(1:N),p(1:S_L),5);
%  figure;
%  [tfd3,orienttfd]=HTFD_new1(extr_Sig2,2,20,84);
% tfsapl(x_o(1:N),tfd3);
% title('Proposed Approach');


% 
%  x=1:N;
%  Sig_r(1:N)=Sig_r(1:N)*WN/2;
%  Sig_r(NA)=x_o(NA);
%  Sig_r=Sig_r(1:N);
%  figure;
%  [tfd3,orienttfd]=HTFD_new1(Sig_r(1:N),2,20,84);
%  title('Stankovic Method');

%  tfsapl(x_o(1:N),tfd3);
% % plot(real(extr_Sig2(1:N)),'r')
% % hold on;
% % plot(real(x_o(1:N)),'b:')
% figure;
% 
% plot(real(extr_Sig2(1:N)),'r')
% hold on;
% plot(real(x_o(1:N)),'b:')
% title('Proposed Approach');
% figure;

plot(real(Sig_r(WN:N)),'r')
hold on;
plot(real(x_o(WN:N)),'b:')
title('Stankovic Method');
display('Stankovic Method')

figure,
plot(real(x_o(WN:N))-real(Sig_r(WN:N)))

% mean(abs(x_o(32:N)-Sig_r(32:N)))
% display('The Proposed')
% mean(abs(x_o(32:N)-extr_Sig2(32:N)))

% I=quadtfd(x_o(1:N),N-1,1,'specx',31,'hamm',N);
% figure
% tfsapl(x_o,abs(I).^2);

% Error in dB in time domain :
Error_DB_time1=10*log10(sum(abs(x_o(WN:N)-Sig_r(WN:N)).^2))
% Error_DB_time2=10*log10(sum(abs(extr_Sig2(WN:N)-Sig_r(WN:N)).^2))


end