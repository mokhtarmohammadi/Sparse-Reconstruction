clear all, close all, clc

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


n=1:1024+256;
N=256;


% Original signal in time domain:
%     x_original=ones(size(n));
   x_original=1.5*exp(j*2*pi*96*n/N+j*24*pi*n.^2/N^2+j*2*pi/6)/N+exp(j*2*pi*24*n/N+j*8*pi*n.^2/N^2+j*2*pi/6)/N;

   
SignalRekon=zeros(1024+N,1); % Initial value of the reconstructed signal in time domain


STFT=[]; STFT_R=[];
it=0;
for tt=0:32:1024     % Signal length (1024+256) with a step of 32
    it=it+1;
    
    N=256;                % Length of the window
    NA=round(3*N/4);     % Number of available samples
    K=48;                % Sparsity
    n=tt+(1:N)';         % Signal index
    
    % Generate the windowed signal (same as x_original, just windowed):
    x=1.5*exp(j*2*pi*96*n/N+j*24*pi*n.^2/N^2+j*2*pi/6).*hanning(N,'periodic')/N;
    x=x+exp(j*2*pi*24*n/N+j*8*pi*n.^2/N^2+j*2*pi/6).*hanning(N,'periodic')/N;

%      x=ones(size(n)).*hanning(N,'periodic');
    
       
    X=fft(x);                    % FT of the original (nonsparse) signal
    STFT=[STFT, X];              % STFT of the original (nonsparse) signal
    

    p=randperm(N); NA_v=p(1:NA); % Random positions of the available samples
    
    
    % Generate the measurement matrix:
    IDFT_mtx=conj(dftmtx(N))/N;     % Full IDFT matrix
    A=IDFT_mtx(NA_v,:);             % Partial IDFT matrix
    
    y=A*X;                          % Available samples
    
 
    % Reconstruction with iterative OMP:
    KB=[];  y0=y; Kstep=1;
    for iter=1:round(K/Kstep)
        
        if iter>1; y=y0-x_R; end            % Remove the reconstructed component from the array
        
        X0=(A*N)'*y;                        % Initial estimate
        [Xv,Kv]=sort(abs(X0));              % Sort the components
        KB=union(KB,Kv(end-Kstep+1:end));   % Take the largest component
        A_K=A(:,KB);                        % Measurement matrix for the largest component
        X_R=zeros(N,1);
        X_R(KB)=pinv(A_K)*y0;               % FT reconstruction of the components
        x_R=A*X_R;                          % Inverse of X_R
    end
    
    STFT_R=[STFT_R, X_R]; % Reconstructed STFT
    
    
    SignalR=ifft(X_R)/16; % Reconstructed window signal
    
    % Reconstructed whole signal (see pages 552 - 554 in DSP book attached)
    SignalRekon(n)=SignalRekon(n)+SignalR;  
    
    
    % Statistical error per one window:
    E_stat(it)=10*log10(sum(abs(X_R(KB)-X(KB)).^2));
    
    KT=1:N; KT(KB)=[]; % Set of non-reconstructed components
    
    C=K*(N-NA)/NA/N; % Value needed for error calculation
    
    % Theoretical error per one window:
    E_theor(it)=10*log10(C*sum(abs(X(KT)).^2));
    
end

% Average of the statistical and theoretical error:
ES=mean(E_stat)
ET=mean(E_theor)

% S-method representation (only used for plotting):
L=15;                       % Length of the window
SM=SM_calc(STFT,L);         % S-method of the original (nonsparse) signal
SM_R=SM_calc(STFT_R,L);     % S-method of the reconstructed signal

% Plotting STFT signals:
figure(1), subplot(211), colormap([0 0 0])
waterfall(abs(STFT).^1'),view([10 50]), axis([0 N 0 32 0 1])
title('Original (nonsparse) STFT')
set(gca,'XTick',[0 100 200],'YTick',[0 16 32])
text(265,5,0,'time'), text(25,-8,0,'frequency')

subplot(212)
waterfall(abs(STFT_R.^1)'),view([10 50]), axis([0 N 0 32 0 1])
title('Reconstructed STFT')
set(gca,'XTick',[0 100 200],'YTick',[0 16 32])
text(265,5,0,'time'), text(25,-8,0,'frequency')


%Plotting S-method signals:
figure(2), subplot(211), colormap([0 0 0])
waterfall(abs(SM).^1'),view([10 50]), axis([0 N 0 32 0 1])
title('Original (nonsparse) SM')
set(gca,'XTick',[0 100 200],'YTick',[0 16 32])
text(265,5,0,'time'), text(25,-8,0,'frequency')

subplot(212)
waterfall(abs(SM_R.^1)'),view([10 50]), axis([0 N 0 32 0 1])
title('Reconstructed SM')
set(gca,'XTick',[0 100 200],'YTick',[0 16 32])
text(265,5,0,'time'), text(25,-8,0,'frequency')


% Comparison of the reconstructed and original signal in the time domain
figure(3),subplot(311)
plot(real(x_original(N:1024)))
title('Original signal in time domain'),xlabel('Time'),ylabel('Amplitude')

subplot(312),plot(real(4*SignalRekon(N:1024)))
title('Reconstructed signal in time domain'),xlabel('Time'),ylabel('Amplitude')

subplot(313),plot(real(x_original(N:1024))'-real(4*SignalRekon(N:1024)))
title('Comparison between the original and the reconstructed signal in time domain')
xlabel('Time'),ylabel('Amplitude')

% Error in dB in time domain :
Error_DB_time=10*log10(sum(abs(x_original(N:1024)'-4*SignalRekon(N:1024)).^2))

