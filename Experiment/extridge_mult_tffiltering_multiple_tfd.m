function [extr_Sig1] = extridge_mult_tffiltering_multiple_tfd(Sig, num, delta,bw,Seg,iii,iter)
th=10;
WIN=length(Sig)/2;
if (isreal(Sig))
    Sig = hilbert(Sig);
end


extr_Sig1=Sig;
extr_Sig1(iii)=0;
mean(abs(Seg-Sig).^2)
%mean(abs(real(Seg-Sig)).^2)
mean(abs(real(Seg)-real(extr_Sig1)).^2)

for ii=1:iter
    Sig=extr_Sig1;
    extr_Sig1(iii)=0;
    
    for i=1:num
        
        [I2,orient]=HTFD_new1(Sig,2,20,104);
        Spec=I2;
        seuilTFD = th*max(Spec(:))/100;
        Spec(Spec<seuilTFD)=0;
        
        c = findridges_new1(Spec,orient,delta);
        
        IF=(c)/(2*length(Spec));
        
        Phase=2*pi*filter(1,[1 -1],IF);
        s_dechirp=exp(-1i*Phase);
        L=5;
        %TF filtering for each sensor
        s1 = Sig.*(s_dechirp);
        s2=fftshift(fft(s1));
        s3=zeros(1,length(Sig));
        s3(WIN-L:WIN+L)=s2(WIN-L:WIN+L);%.*hamming(2*L+1)';
        s2(WIN-L:WIN+L)=0;
        extr_Sig=ifft(ifftshift(s3)).*conj(s_dechirp);
        s2=ifft(ifftshift(s2)).*conj(s_dechirp);
        Sig=s2;%-extr_Sig(iii);
        extr_Sig1(iii)=extr_Sig1(iii)+extr_Sig(iii);
    end
    mean(abs(extr_Sig1-Seg).^2)
    %  mean(abs(real(extr_Sig1)-real(Seg)).^2)
    
end

% save extridge Spec1 Spec2 Spec3 Spec4 Sig1 Sig2 Sig3 Sig4;