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
addpath([upperPath '\TFSA7'])
addpath([upperPath '\Core\AOK'])
addpath([upperPath '\TFTB2'])

%% Signal Parameters and Loading
 tr = 5;
SampFreq = 128;
fs=SampFreq;
t = 0:1/SampFreq:1-1/SampFreq;
N=SampFreq;
n = 0:127;
% M=64;
%M=32;
iter=10;
Sig1 = exp(1i*(2*pi*(0.1*n.^3)/(128*128))+1i*(2*pi*(0.05*n))+1i*(2*pi*(0.05*n.^2/128))); %300t»???150t
Sig2= exp(1i*(2*pi*(0.1*n.^3)/(128*128))+1i*(2*pi*(0.15*n))+1i*(2*pi*(0.05*n.^2/128)));
%  Sig2= exp(1i*(-2*pi*(0.1*n.^3)/(128*128))+1i*(2*pi*(0.45*n))-1i*(2*pi*(0.05*n.^2/128)));
IF_O(:,1)=(.3/(128*128))*n.^2+.05+(.1/128)*n;
IF_O(:,2)=(.3/(128*128))*n.^2+.15+(.1/128)*n;

Sig4 =exp(1i*(2*pi*(50*t - 5*t.^3)));
 %Sig1=exp(1i*(2*pi*(15*t)));
 %Sig4=exp(1i*(2*pi*(30*t)));
% Sig=Sig1.*hamming(128)'+Sig2.*(64-abs([-64:1:63]))/64;
%Sig=real(Sig);
%Sig=Sig1+Sig4;
Sig=Sig1+Sig2;

seg=Sig;
%Sig(80:90)=0;
%Sig(50:65)=0;

%iiii=1;
 %Sig=y;
%  load Data1;
%%%%%%%%%%%STFT%%%%%%%%%%%%%%%%%%%%%
window = 128;
Nfrebin = length(Sig);
aa=1;
f = linspace(0,SampFreq/2,length(Sig));
%% ridge extraction and fitting
num = 2; % the number of the components
delta = 20/4;
N_S=2;
ii=1; 
M=[16 32 48 64];
for jjjj=1:4
for K=1:N_S
    [x Sig iii iiii]=signal_2(SampFreq,M(jjjj));
       %         [extr_Sig1,fidexmult,tfdv] =extridge_mult_tfds(Sig,'ADTFD', SampFreq, num, delta, orderIF1,bw1,Nfrebin,window,alpha,iiii);
[extr_Sig1,fidexmult] =extridge_mult_tffiltering_multiple_tfd_viterbi(Sig, num, delta,3,x,iii,3);
% extr_Sig2=extr_Sig1*mean(abs(Sig(iiii)))/mean(abs(extr_Sig1(iiii)));
% extr_Sig2=extr_Sig2*exp(1i*angle(sum(extr_Sig2(iiii).*conj(Sig(iiii)))));

adtfd_precision(K)=10*log10(sum(abs(extr_Sig1-x).^2)/length(x));
% adtfd_precision
%%----------------------------------

Precision=1000; % Precision in dB
yr1=GradRec(Sig,iii,{@MySTFT,150},delta,60);
grad_precision(K)=10*log10(sum(abs(yr1-x).^2)/length(x));
grad_precision
%%----------------------
yr2=MyOMP(Sig,iii,N);
precision_OMP(K)=10*log10(sum(abs(yr2-x).^2)/length(x));
precision_OMP
%%---------------
% yr=MyBP(Sig,N);
% yr3=MyIterBP(Sig,N);
% precision_BP(K)=10*log10(sum(abs(yr3-x).^2)/length(x));
% precision_BP
end
adtfd_precision_mean(jjjj)=mean(adtfd_precision);
grad_precision_mean(jjjj)=mean(grad_precision);
OMP_precision_mean(jjjj)=mean(precision_OMP);
% BP_precision_mean(jjjj)=mean(precision_BP);
end

%end
%%----------------------------------------
% yr=CTD1(x,iiii,N-M);
%for K=1:N_S
       % [x Sig iii iiii Nfrebin]=signal_1(SampFreq,N,M);
%yr=MyOMP(iii,N);
%Obtained_precision_OMP=10*log10(sum(abs(yr-x).^2)/length(x))
%end
%mean(Obtained_precision_OMP)

%%------------------------------------
%for K=1:N_S
        %[x Sig iii iiii Nfrebin]=signal_1(SampFreq,N,M)

%yr=MyBP(iii,N);
%Obtained_precision_BP=10*log10(sum(abs(yr-x).^2)/length(x))
%mean(Obtained_precision_BP)

% NS=32;
figure;
    plot(M,adtfd_precision_mean,'-co','linewidth',4);
    hold on;
    plot(M,grad_precision_mean,':gs','linewidth',4);
   hold on;
     plot(M,OMP_precision_mean,'-.b+','linewidth',4);
    hold on;
%     plot(M,BP_precision_mean,'--md','linewidth',4);
    legend(' ADTFD','Gradient based',' Extended OMP',' Iterative Thresholding')
    xlabel('Number of missing samples ');
    ylabel('Mean Square Error (dB)');
    