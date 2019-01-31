%%  Author:
%     Mokhtar Mohammadi
%In this code we assume that the user add in this path the TFSAP toolboox.
%% Locating and Adding Functions Directory
clear all;
close all;
currentDirectory = pwd;
[upperPath, ~, ~] = fileparts(currentDirectory);
addpath([upperPath '\Core'])
addpath([upperPath '\TFSA7'])
addpath([upperPath '\Core\AOK'])
addpath([upperPath '\TFTB2'])

%% Signal Parameters and Loading
tr=1;
SampFreq = 142000;
fs=SampFreq;
% t = 0:1/SampFreq:1-1/SampFreq;
 N=400;
% n = 0:N-1;
M=100;
%% Loading
load bat1;
fs = 142000; tr = 5;
Sig=(bat1(:)');
seg=Sig;
p = randperm(N);
y = Sig;
x=y;
y(p(1:M)) = 0;
y0=y;               % initial y
iii=p(1:M);
count=1;
for jj=1:400
   if isempty(find(iii==jj, 1))
       iiii(count)=jj;
       count=count+1;
   end
end
%iiii=1;
 Sig=hilbert(y);
 %Sig=y;
%  load Data1;
%%%%%%%%%%%STFT%%%%%%%%%%%%%%%%%%%%%
Ratio = 0;
window = 128;
Nfrebin = length(Sig);
aa=1;
f = linspace(0,SampFreq/2,length(Sig));
%% ridge extraction and fitting
bw1 = SampFreq/60;%
orderIF1 = 50; % use the Fourier model to smooth the detected ridge curves£»orderIF1 could be larger than the following orderIF
num = 3; % the number of the components
delta = 20/4;
alpha = 5;
[tfd1,orienttfd]=HTFD_new1(seg,2,20,96);
%  [tfd,Orient]=HTFD_new_post(tfd,3,8,84);
%figure; SetFigDef(16,9,'Times',20);tfsapl(x,tfd1(:,1:1:end),'SampleFreq',fs,'Title','(a)','TFfontSize' , 20,'plotfn','pcolor');
figure; imagesc(tfd1);
[tfd2,orienttfd]=HTFD_new1(Sig,2,20,96);
%figure; SetFigDef(16,9,'Times',20);tfsapl(x,tfd2(:,1:1:end),'SampleFreq',fs,'Title','(b)','TFfontSize' , 20,'plotfn','pcolor');
figure;imagesc(tfd2);
%[extr_Sig1,fidexmult] = extridge_mult_new_modifiedTFFILTER(Sig,tfd,orienttfd, num, delta,iii);

%[extr_Sig1,fidexmult,tfdv] =extridge_mult_new(Sig, SampFreq, num, delta, orderIF1,bw1,Nfrebin,window,alpha,iiii);
[extr_Sig,fidexmult,tfdv] =extridge_mult_tffiltering(Sig, SampFreq, num, delta,bw1,iiii,iii);
extr_Sig2=extr_Sig*mean(abs(Sig(iiii)))/mean(abs(extr_Sig(iiii)));
%extr_Sig1=extr_Sig1*exp(1i*angle(sum(extr_Sig1.*conj(Sig))));
extr_Sig2=extr_Sig2*exp(1i*angle(sum(extr_Sig2(iiii).*conj(Sig(iiii)))));

[tfd3,orienttfd]=HTFD_new1(extr_Sig2,2,20,96);
% [tfdr,orienttfd]=F_ADTFD(extr_Sig2,2,30,80,[0 15 175 180]);
%figure; SetFigDef(16,9,'Times',20);tfsapl(extr_Sig2,tfd3(:,1:1:end),'SampleFreq',fs,'Title','(c)','TFfontSize' , 20,'plotfn','pcolor');
figure; imagesc(tfd3);
%  [tfdr,Orient]=HTFD_new_post(tfdr,3,8,84);
% save Data1 Sig x;