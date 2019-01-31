clear all, close all, clc
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
WN=64;
n=0:127+WN;
Sig1 = exp(1i*(2*pi*(0.1*n.^3)/(128*128))+1i*(2*pi*(0.05*n))+1i*(2*pi*(0.05*n.^2/128))); %300t»???150t
Sig2= exp(1i*(2*pi*(0.051*n.^3)/(128*128))+1i*(2*pi*(0.2*n))+1i*(2*pi*(0.05*n.^2/128)));
Sig2= exp(1i*(-2*pi*(0.1*n.^3)/(128*128))+1i*(2*pi*(0.45*n))-1i*(2*pi*(0.05*n.^2/128)));

%Sig1 = exp(1i*(2*pi*(0.1*n))); %300t»???150t
%Sig2= exp(1i*(2*pi*(0.3*n)));
n=0:127;



x_o=Sig1+Sig2;

n=0:255;
for jj=1:WN
    IDFT_mtx(jj,:)=exp(-1i*2*pi*(jj-1)*n/256);     % Full IDFT matrix
end
x_p=x_o;
IDFT_mtx=conj(IDFT_mtx)/WN;
%x_p=[zeros(1,WN/2) x_o zeros(1,WN/2)];
x_s=x_p;

x_p=[zeros(1,WN/2) x_o ];
p = randperm(length(x_p));
S_L=56;
 p(1:S_L)=[10:30 60:80  90:103 ];
 
%S_L=84;% Sparsity level

x_s(p(1:S_L))=0;
x_s=[zeros(1,WN/2) x_s ];
NA=find(x_s~=0);

STFT=[]; STFT_R=[];
it=0;

for tt=0:1:length(x_p)-WN-WN/2 -1    % Signal length (1024+256) with a step of 32
    it=it+1;
    
    
    K=2;                % Sparsity
    %n=tt+(1:N)';         % Signal index
    
    % Generate the windowed signal:
    %     x=1.5*exp(1i*2*pi*96*n/N+1i*24*pi*n.^2/N^2+1i*2*pi*rand).*blackmanharris(N,'periodic')/N;
    %     x=x+exp(1i*2*pi*24*n/N+1i*8*pi*n.^2/N^2+1i*2*pi*rand).*blackmanharris(N,'periodic')/N;
    %     n=n(1:32);
    %     xx=1.5*exp(1i*2*pi*96*n/N+1i*24*pi*n.^2/N^2+1i*2*pi*rand)+exp(1i*2*pi*24*n/N+1i*8*pi*n.^2/N^2+1i*2*pi*rand);
    %     x_s=[x_s; xx];
    %
    x_w=x_p(tt+1:tt+WN).'.*(blackmanharris(WN))/WN;
    
    
    X=(fft(x_w,256));                    % FT of the original (nonsparse) signal
   % X=X(1:128);
    %X=(fft(x_w,128));                    % FT of the original (nonsparse) signal
    
    STFT=[STFT, X(1:128)];              % STFT of the original (nonsparse) signal
    
    
    NA_v=NA(and(NA<tt+WN,NA>tt))-tt; % Random positions of the available samples
    
    
    % Generate the measurement matrix:
    %IDFT_mtx=conj(dftmtx(length(x_w)))/length(x_w);     % Full IDFT matrix
    A=IDFT_mtx(NA_v,:);             % Partial IDFT matrix
    
    
    y=A*X;                          % Available samples
    
    
    
    % Reconstruction with iterative OMP:
    KB=[];  y0=y; Kstep=1;
    for iter=1:round(K/Kstep)
        
        if iter>1; y=y0-x_R; end            % Remove the reconstructed component from the array
        
        X0=(A*WN)'*y;                        % Initial estimate
        [Xv,Kv]=sort(abs(X0));              % Sort the components
        KB=[KB; Kv(end-Kstep+1:end)];       % Take the largest component
        A_K=A(:,KB);                        % Measurement matrix for the largest component
        X_R=zeros(256,1);
        X_R(KB)=pinv(A_K)*y0;               % FT reconstruction of the components
        x_R=A*X_R;                          % Inverse of X_R
    end
    
    STFT_R=[STFT_R, X_R(1:128)]; % Reconstructed STFT
    
    % Statistical error per one window:
    %     E_stat(it)=10*log10(sum(abs(X_R(KB)-X(KB)).^2));
    %
    %     KT=1:length(x_w); KT(KB)=[]; % Set of non-reconstructed components
    %
    %     C=K*(WN-NA)/NA/WN; % Value needed for error calculation
    %
    %     % Theoretical error per one window:
    %     E_theor(it)=10*log10(C*sum(abs(X(KT)).^2));
    %
end
figure(3);
% sig_r=istft(STFT_R,blackmanharris(N,'periodic')/N, hamming(N,'periodic')/N, 32, 256, 2);
% plot(real(sig_r)*4,'r')
% hold on; plot(real(x_s(WN/2+1:128+WN/2-1+1)),'b:');
% Average of the statistical and theoretical error:
%ES=mean(E_stat);
%ET=mean(E_theor);

% S-method representation (only used for plotting):

% Plotting STFT signals:
%tfsapl(x_o,abs(STFT).^2);
%figure
%tfsapl(x_o,abs(STFT_R).^2);
% [extr_Sig2]=extridge_mult_tffiltering_multiple_tfd_viterbi(x_s(WN/2+1:128+WN/2-1+1), 2, 5,3,x_o(1:128),p(1:S_L),5);
% figure;
% [tfd3,orienttfd]=HTFD_new1(extr_Sig2,2,20,84);
% tfsapl(x_o,tfd3);
figure
tfsapl(x_o,abs(STFT_R).^2);
I=quadtfd(x_o(1:128),128-1,1,'specx',31,'hamm',128);
figure
tfsapl(x_o,abs(I).^2);

%
% %Plotting S-method signals:
% figure(2), subplot(211), colormap([0 0 0])
% waterfall(abs(SM).^1'),view([10 50]), axis([0 N 0 32 0 1])
% title('Original (nonsparse) SM')
% set(gca,'XTick',[0 100 200],'YTick',[0 16 32])
% text(265,5,0,'time'), text(25,-8,0,'frequency')
%
% subplot(212)
% waterfall(abs(SM_R.^1)'),view([10 50]), axis([0 N 0 32 0 1])
% title('Reconstructed SM')
% set(gca,'XTick',[0 100 200],'YTick',[0 16 32])
% text(265,5,0,'time'), text(25,-8,0,'frequency')
%
%
