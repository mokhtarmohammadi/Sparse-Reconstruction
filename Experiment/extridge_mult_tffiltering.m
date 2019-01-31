function [extr_Sig1,fidexmult, tfdv] = extridge_mult_tffiltering(Sig, SampFreq, num, delta,bw,iiii,iii)
fs=SampFreq;
th=10;
WIN=length(Sig)/2;
if (isreal(Sig))
Sig = hilbert(Sig);
end

fidexmult = zeros(num,length(Sig));
tfdv = zeros(num,length(Sig));

extr_Sig1=Sig;
for i = 1:num
f = linspace(0,SampFreq/2,length(Sig));
%I2=HTFD_new1(Sig,2,20,128);
Sig1=Sig;
%Sig1(iii)=0;
I2=HTFD_new1(Sig1,2,15,64);
I2(I2<0)=0;
Spec=I2;
seuilTFD = th*max(Spec(:))/100;
    Spec(Spec<seuilTFD)=0;

c = findridges(Spec,delta);
            
                 IF=(c)/(2*length(Spec));
                Phase=2*pi*filter(1,[1 -1],IF);
                s_dechirp=exp(-1i*Phase);
                L=3;
                %TF filtering for each sensor
                s1 = Sig.*(s_dechirp);
                s2=fftshift(fft(s1));
                s3=zeros(1,length(Sig));
                 s3(WIN-L:WIN+L)=s2(WIN-L:WIN+L);
                s2(WIN-L:WIN+L)=0;
                extr_Sig=ifft(ifftshift(s3)).*conj(s_dechirp);
                s2=ifft(ifftshift(s2)).*conj(s_dechirp);
                Sig=s2;%-extr_Sig(iii);
               extr_Sig1(iii)=extr_Sig1(iii)+extr_Sig(iii);
 % remove the extracted signal component so that other ridge curves with smaller energies can be extracted in the subsequent iterations
end
% save extridge Spec1 Spec2 Spec3 Spec4 Sig1 Sig2 Sig3 Sig4;