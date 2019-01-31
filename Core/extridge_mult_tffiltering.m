function [extr_Sig1,fidexmult, tfdv] = extridge_mult_tffiltering(Sig, SampFreq, num, delta,bw,iiii,iii)
fs=SampFreq;
th=10;
WIN=length(Sig)/2;
if (isreal(Sig))
Sig = hilbert(Sig);
end

orderamp = round(bw*length(Sig)/SampFreq);%Fourier order for characterizing signal amplitude
fidexmult = zeros(num,length(Sig));
tfdv = zeros(num,length(Sig));
extr_Sig1=zeros(1,length(Sig));
extr_Sig1=Sig;
for i = 1:num
f = linspace(0,SampFreq/2,length(Sig));
I2=HTFD_new1(Sig,2,20,128);
% I2=HTFD_new1(Sig,3,8,64);
I2(I2<0)=0;
Spec=I2;
seuilTFD = th*max(Spec(:))/100;
    Spec(Spec<seuilTFD)=0;
% if i==1
% Spec1=Spec;
% Sig1=Sig;
% % figure; SetFigDef(16,9,'Times',20);tfsapl(Sig,Spec(:,1:1:end),'SampleFreq',fs,'Title','(b)','TFfontSize' , 20,'plotfn','pcolor');
% elseif i==2
% Spec2=Spec;
% Sig2=Sig;
% % figure; SetFigDef(16,9,'Times',20);tfsapl(Sig,Spec(:,1:1:end),'SampleFreq',fs,'Title','(c)','TFfontSize' , 20,'plotfn','pcolor');
% elseif i==3 
% Spec3=Spec;
% Sig3=Sig;
% %     figure; SetFigDef(16,9,'Times',20);tfsapl(Sig,Spec(:,1:1:end),'SampleFreq',fs,'Title','(d)','TFfontSize' , 20,'plotfn','pcolor');
% else
%     th=50;
%     seuilTFD = th*max(Spec(:))/100;
%     Spec(Spec<seuilTFD)=0;
% %     figure; SetFigDef(16,9,'Times',20);tfsapl(Sig,Spec(:,1:1:end),'SampleFreq',fs,'Title','(e)','TFfontSize' , 20,'plotfn','pcolor');
% Sig4=Sig;
% Spec4=Spec;
% end

c = findridges(Spec,delta);
            
                 IF=(length(Sig)-c)/(2*length(Sig));
                 IF=(c)/(2*length(Sig));
                Phase=2*pi*filter(1,[1 -1],IF);
                s_dechirp=exp(-1i*Phase);
                L=4;
                %TF filtering for each sensor
                s1 = Sig.*(s_dechirp);
                s2=fftshift(fft(s1));
                s3=zeros(1,length(Sig));
%                 s3(128-L:128+L)=s2(128-L:128+L);
%                 s2(128-L:128+L)=0;
                 s3(WIN-L:WIN+L)=s2(WIN-L:WIN+L);
                s2(WIN-L:WIN+L)=0;
                extr_Sig=ifft(ifftshift(s3)).*conj(s_dechirp);
                s2=ifft(ifftshift(s2)).*conj(s_dechirp);
                Sig=s2;%-extr_Sig(iii);
               extr_Sig1(iii)=extr_Sig1(iii)+extr_Sig(iii);
 % remove the extracted signal component so that other ridge curves with smaller energies can be extracted in the subsequent iterations
end
% save extridge Spec1 Spec2 Spec3 Spec4 Sig1 Sig2 Sig3 Sig4;