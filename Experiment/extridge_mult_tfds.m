function [extr_Sig1,fidexmult, tfdv] = extridge_mult_tfds(Sig,TFD, SampFreq, num, delta, orderIF,bw,Nfrebin,window,alpha,iii)
% Extract ridges for multi-component signals.
% In each iteration,the signal component associated with the extrated ridge is
% reconstructed by the ICCD and then removed from the original signal so
% that the ridge curves of other signal components with smaller energies
% can be extracted in the subsequent iterations.
%%%%%%%%%%%%%%%%%%%%%%%    input      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sig��measured signal,a row vector
% SampFreq: sampling frequency
% num: the number of the signal components
% delta��maximum allowable frequency variation between two consecutive points
% orderIF: the order of the Fourier model used for smoothing the extracted ridge curves
% bw��the bandwidth of the ICCD (unit��Hz); herein the ICCD can be regarded as a time-frequency filtering technique
% Nfrebin,window are two parameters for implementing the STFT
% alpha��Tikhonov regularization parameter for ICCD.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  output   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fidexmult: obtained ridge index for multi-component signals
% tfdv: the corresponding ridge magnitude 
fs=SampFreq;
if (isreal(Sig))
Sig = hilbert(Sig);
end

orderamp = round(bw*length(Sig)/SampFreq);%Fourier order for characterizing signal amplitude
fidexmult = zeros(num,length(Sig));
tfdv = zeros(num,length(Sig));
extr_Sig1=zeros(1,length(Sig));
for i = 1:num
%[Spec,f] = STFT(Sig(:),SampFreq,Nfrebin,window); %STFT
%f = linspace(-SampFreq/2,SampFreq/2,nLevel);
f = linspace(0,SampFreq/2,length(Sig));
switch TFD
    case 'WVD'
        I2 = quadtfd( Sig, length(Sig)-1, 1, 'wvd');
    case 'SPEC'
        I2= quadtfd( Sig, length(Sig)-1, 1, 'specx', 23, 'hann',128);
    case 'EMBD'
        I2= quadtfd(Sig,length(Sig)-1, 1, 'emb', 0.1, 0.1);
    case 'CKD'
        I2= cmpt( Sig, 'ckd', 1, 0.1, 0.1);
    case 'SM'
        I2= Adaptive_S_method(Sig,23,'hann');
    case 'AOK'
        I_max_new=adaptive_optimal_tfd(Sig);
        I2 = real(I_max_new);
    case 'RGK'
        I2= rgk(Sig,2);
    case 'ADTFD'
        I2=HTFD_new1(Sig,2,30,80);
end
%I2=Sparse_adtfd(Sig,2,24,64*2,num-i+1);
%I2=tfr_stft_high(Sig);
%I2 = quadtfd(Sig, length(Sig)-1, 1, 'specx',95,'hamm',length(Sig));
I2(I2<0)=0;
%  I2=imresize(I2,[length(Sig)/2 length(Sig)]);
% I3=zeros(length(Sig),length(Sig));
% I3(end/2+1:end,:)=I2;
% 
% Spec=I3;

Spec=I2;

% imagesc(Spec)
c = findridges(Spec,delta);
[extr_Sig,~,IFfit] = ICCD(Sig,SampFreq,f(c),orderIF,orderamp,alpha);
extr_Sig1=extr_Sig1+extr_Sig;
findex = zeros(1,length(Sig));
for j = 1:length(Sig)
   [~,findex(j)] = min(abs(f - IFfit(1,j)));
   tfdv(i,j) = abs(Spec(findex(j),j));
end
fidexmult(i,:) = findex;
Sig(iii) = Sig(iii) - extr_Sig(iii); % remove the extracted signal component so that other ridge curves with smaller energies can be extracted in the subsequent iterations
end
