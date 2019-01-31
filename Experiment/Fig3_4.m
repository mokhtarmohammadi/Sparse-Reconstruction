%%  Author:
%     Mokhtar Mohammadi
%In this code we assume that the user add in this path the TFSAP toolboox.
%% Locating and Adding Functions Directory
clear all;
close all;
addpath('D:\tfsa_5-5\windows\win64_bin');
currentDirectory = pwd;
[upperPath, ~, ~] = fileparts(currentDirectory);
addpath([upperPath '\Core'])
% addpath([upperPath '\TFSA7'])
% addpath([upperPath '\Core\AOK'])
% addpath([upperPath '\TFTB2'])

%% Signal Parameters and Loading
 tr = 5;
SampFreq = 128;
fs=SampFreq;
t = 0:1/SampFreq:1-1/SampFreq;
N=SampFreq;
M=64;
WN=64; N=128;
n=0:N-1+WN*0;

%M=32;
iter=10;
Sig1 = exp(1i*(2*pi*(0.1*n.^3)/(128*128))+1i*(2*pi*(0.05*n))+1i*(2*pi*(0.05*n.^2/128))); %300t»???150t
Sig2= exp(1i*(2*pi*(0.1*n.^3)/(128*128))+1i*(2*pi*(0.15*n))+1i*(2*pi*(0.05*n.^2/128)));
%  Sig2= exp(1i*(-2*pi*(0.1*n.^3)/(128*128))+1i*(2*pi*(0.45*n))-1i*(2*pi*(0.05*n.^2/128)));
IF_O(:,1)=(.3/(128*128))*n.^2+.05+(.1/128)*n;
IF_O(:,2)=(.3/(128*128))*n.^2+.15+(.1/128)*n;
x_o=Sig1+Sig2;
x_p=[zeros(1,WN) x_o zeros(1,WN)];
x_s=x_p;
S_L=56;
%p = randperm(N);
% load Data2;
p(1:S_L)=[10:30 60:80  90:103 ];

% M=48;
% p(1:M)=[11:30 60:75  90:101 ];%84:84+7 100:107];
% 
%  M=56;
%  p(1:M)=[10:30 60:80  90:103 ];

%  M=62;
%  p(1:M)=[10:32 60:82  90:105 ];
 
 % S_L=64;
 % p(1:M)=[9:32 60:83  90:105 ];

%  M=64;
 % p(1:M)=[10:30 60:80  90:111 ];

% M=64;
% p(1:M)=[10:36 58:79  95:109 ];
%p(1:M)=[10:17 24:31 40:47 60:67 80:87 100:107];%84:84+7 100:107];
 %M=32;
 %p(1:M)=[10:25 80:95];%84:84+7 100:107];
x_s(p(1:S_L)+WN)=0;
% x_s=[ x_s ];

NA=find(x_s~=0);

pp=find(x_s==0);
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
% M4=mean(abs((Sig_r(WN+1:N+WN))-(x_p(WN+1:N+WN))).^2)

%%%%%%%%%%%STFT%%%%%%%%%%%%%%%%%%%%%
 Sig=x_s(WN+1:N+WN);
% window = 128;
% Nfrebin = length(Sig);
% aa=1;
% f = linspace(0,SampFreq/2,length(Sig));
% %% ridge extraction and fitting
% num = 2; % the number of the components
% delta = 20/4;
 [tfd1,~]=HTFD_new1(Sig_r(WN+1:N+WN),2,50,84);
% % figure; tfsapl(real(Sig),tfd1);
%%%%%%%%%%%%%
[tfd,~]=HTFD_new1(x_s(WN+1:N+WN),2,50,84);
figure;SetFigDef(16,9,'Times',20); tfsapl(real(x_s(WN+1:N+WN)),tfd,'YLabel','Sample Number','Title','(a)', 'TFfontSize' , 20,'grayscale','on');
IF_est= ADTFD_IF_estimation_viterbi_modified(x_s(WN+1:N+WN), 2,5);
t=0:127;
figure;
SetFigDef(16,9,'Times',20);	    
% set(gcf,'Color','w'); 
plot(IF_est(:,1:4:end),t(1:4:end),'o',IF_O',t,'-','linewidth',2);
ylabel('Sample Number');
xlabel('Instantaneous Frequency (Hz)');
title('(b)');
axis([0 0.5 0 128]);
% M0=mean(abs(extr_Sig1-Seg).^2)
% ans=  0.9460  0.6657   0.1386
      
%[extr_Sig2]=extridge_mult_tffiltering_new(Sig, num, delta,3,seg,iii,iter);
 %  display('second method');
%[extr_Sig2]=extridge_mult_tffiltering_multiple_tfd(Sig, num, delta,3,x_p(WN+1:N+WN),p(1:S_L),3);
%  [extr_Sig2 IF_est]=extridge_mult_tffiltering_multiple_tfd_viterbi_visual(Sig, num, delta,3,seg,iii,5,IF_O);
  %[extr_Sig2]=extridge_mult_tffiltering_multiple_tfd_viterbi(x_s(WN+1:N+WN), 2, 5,8,x_p(WN+1:N+WN),p(1:S_L),5);

[extr_Sig2  IF_est]=extridge_mult_tffiltering_multiple_tfd_viterbi_visual(x_s(WN+1:N+WN), 2, 5,5,x_p(WN+1:N+WN),p(1:S_L),5,IF_O);
  
[tfd3,orienttfd]=HTFD_new1(extr_Sig2,2,50,84);
figure;SetFigDef(16,9,'Times',20); tfsapl(real(extr_Sig2),tfd3,'YLabel','Sample Number','Title','(e)', 'TFfontSize' , 20, 'grayscale','on');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=0:127;
figure
SetFigDef(16,9,'Times',20);	    
% set(gcf,'Color','w'); 
plot(IF_est(:,1:4:end),t(1:4:end),'o',IF_O',t,'-','linewidth',2);
ylabel('Sample Number');
xlabel('Instantaneous Frequency (Hz)');
title('(f)');
axis([0 0.5 0 128]);
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%
% figure;SetFigDef(16,9,'Times',20); tfsapl(Sig_r(WN+1:N+WN),tfd1,'YLabel','Sample Number','Title','(a)', 'TFfontSize' , 20,'grayscale','on');

% sparsity_aware(x_o,p(1:M),M);
%%%%%%%%%%%%%%%%%
Precision=1000; % Precision in dB
yr1=GradRec(Sig,p(1:S_L),Precision,max(abs(Sig)));
% yr2=MyOMP(Sig,iii,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;SetFigDef(16,9,'Times',20);
plot(real(Sig),'r:','linewidth',2);
hold on; 
% figure;SetFigDef(16,9,'Times',20);
plot(real(x_p(WN+1:N+WN)),'b');
% subplot(5,1,3)
ylabel('Amplitude');
xlabel('Sample Number');

title('(a)');
axis([0 128 -2 2]);
%%%%%%%%%%%%%%%%%%
figure;SetFigDef(16,9,'Times',20);
% subplot(5,1,5)
plot(real(extr_Sig2),'r:','linewidth',2);
hold on; 
% figure;SetFigDef(16,9,'Times',20);
plot(real(x_p(WN+1:N+WN)),'b');
ylabel('Amplitude');
xlabel('Sample Number');

title('(b)');
axis([0 128 -2 2]);


% subplot(5,1,4)
figure;SetFigDef(16,9,'Times',20);
plot(real(yr1),'r:','linewidth',2)
hold on;
plot(real(x_p(WN+1:N+WN)),'b');
ylabel('Amplitude');
xlabel('Sample Number');

title('(c)');
axis([0 128 -2 2]);
% subplot(5,1,2)
figure;SetFigDef(16,9,'Times',20);
plot(real(Sig_r(WN+1:N+WN)),'r:','linewidth',2);
hold on;
plot(real(x_p(WN+1:N+WN)),'b');
ylabel('Amplitude');
xlabel('Sample Number');
title('(d)');
axis([0 128 -2 2]);
% title('(g)');
M2=mean(abs((extr_Sig2)-(x_p(WN+1:N+WN))).^2)
M3=mean(abs((yr1)-(x_p(WN+1:N+WN))).^2)
M4=mean(abs((Sig_r(WN+1:N+WN))-(x_p(WN+1:N+WN))).^2)
%%%%%%%%%%%%%%%
% M1 =
% 
%     0.4009
% M2 =
% 
%     0.0636
% M3 =
% 
%     0.9941
% 
% M4 =
%     1.0773

