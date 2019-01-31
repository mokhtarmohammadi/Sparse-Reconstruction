function x=MyBPold(p,N)
%     N=128; % window width
    win=[0 (hanning(N-1))'];   Ew=sum(win.^2)/N;    
    STFT=[];
    for n=0:1:255
      x=signal_2(n,N);
      x(p(:))=0;
      STFT=[STFT,fftshift(fft(fftshift([zeros(1,N/2), x.*win, zeros(1,N/2)]))).'];
    end
    STFT=STFT(N+1:2*N,:); %STFT positive frequencies
    STFT=STFT/sqrt(Ew)/N; %Normalization for the plot

%Constant L
% L=10; %Correction terms: spectrogram L=0, WD L=length(x)/2-1, L=10 Suggested
% STFT=SM_calc(STFT,10);
M=256;
n=0:M-1;
for ii=0:M-1
    A1(ii+1,:)=exp(1i*ii*2*pi*n/M)/M;
end
for n=1:1:255
    y=STFT(:,n);
% Initial transform vector
x0=A1'*y;

%Running the recovery Algorithm from L1-magic
% xp=l1eq_pd(x0,A1,[],y,1e-3);
% xp= perform_l1_recovery(A1,y,'omp');
 xp= perform_l1_recovery(A1,y,'bp');
%  xp=real(OMP(A1,y,8));
%recovered signal in time domain
xprec(n,:)=(ifft(fftshift(xp)));
end
x=sum(xprec)/N;
