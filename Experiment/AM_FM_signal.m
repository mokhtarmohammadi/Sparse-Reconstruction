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
M=64;
%M=32;
iter=10;
Sig1 = exp(1i*(2*pi*(0.1*n.^3)/(128*128))+1i*(2*pi*(0.05*n))+1i*(2*pi*(0.05*n.^2/128))); %300t»???150t
Sig2= exp(1i*(2*pi*(0.1*n.^3)/(128*128))+1i*(2*pi*(0.15*n))+1i*(2*pi*(0.05*n.^2/128)));
% Sig2= exp(1i*(-2*pi*(0.1*n.^3)/(128*128))+1i*(2*pi*(0.45*n))-1i*(2*pi*(0.05*n.^2/128)));


 Sig4 =exp(1i*(2*pi*(50*t - 5*t.^3)));
 %Sig1=exp(1i*(2*pi*(15*t)));
 %Sig4=exp(1i*(2*pi*(30*t)));
Sig=Sig1.*hamming(128)'+Sig2.*(64-abs([-64:1:63]))/64;
%Sig=real(Sig);
%Sig=Sig1+Sig4;
Sig=Sig1+Sig2;
seg=Sig;
%Sig(80:90)=0;
%Sig(50:65)=0;

p = randperm(N);
M=48;
p(1:M)=[11:30 60:75  90:101 ];%84:84+7 100:107];

 M=56;
 p(1:M)=[10:30 60:80  90:103 ];

%  M=62;
%  p(1:M)=[10:32 60:82  90:105 ];
 
 % M=64;
%  p(1:M)=[9:32 60:83  90:105 ];

%  M=64;
 % p(1:M)=[10:30 60:80  90:111 ];

% M=64;
% p(1:M)=[10:36 58:79  95:109 ];
%p(1:M)=[10:17 24:31 40:47 60:67 80:87 100:107];%84:84+7 100:107];
 %M=32;
 %p(1:M)=[10:25 80:95];%84:84+7 100:107];

y = Sig;
x=y;
y(p(1:M)) = 0;
y0=y;               % initial y
iii=p(1:M);

count=1;
for jj=1:N
   if isempty(find(iii==jj, 1))
       iiii(count)=jj;
       count=count+1;
   end
end
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
Sig=y;
[tfd1,~]=HTFD_new1(Sig,2,20,84);
% figure; tfsapl(real(Sig),tfd1);
figure;SetFigDef(16,9,'Times',20); tfsapl(real(Sig),tfd1,'Title','(a)', 'TFfontSize' , 20);

%[extr_Sig2]=extridge_mult_tffiltering_new(Sig, num, delta,3,seg,iii,iter);
%[extr_Sig2]=extridge_mult_tffiltering_viterbi(Sig, num,3, seg,iii,iter);
 %  display('second method');
%[extr_Sig2]=extridge_mult_tffiltering_multiple_tfd(Sig, num, delta,3,seg,iii,10);
 [extr_Sig2]=extridge_mult_tffiltering_multiple_tfd_viterbi(Sig, num, delta,3,seg,iii,5);

[tfd3,orienttfd]=HTFD_new1(extr_Sig2,2,20,84);

figure; tfsapl(real(extr_Sig2),tfd3);


figure;
subplot(2,1,2)
plot(real(extr_Sig2),'r');
hold on; plot(real(seg),'b:')
subplot(2,1,1)
plot(real(Sig),'r');
hold on; plot(real(seg),'b:')
mean(abs((extr_Sig2)-(seg)).^2)
%  [tfdr,Orient]=HTFD_new_post(tfdr,3,8,84);
% save Data1 Sig x;