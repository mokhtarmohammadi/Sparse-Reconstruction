function [extr_Sig1 IF_est] = extridge_mult_tffiltering_multiple_tfd_viterbi_visual(Sig, num, delta,bw,Seg,iii,iter,IF_O)
WIN=length(Sig)/2;
if (isreal(Sig))
    Sig = hilbert(Sig);
end


extr_Sig1=Sig;
extr_Sig1(iii)=0;
%mean(abs(real(Seg-Sig)).^2)
Sigg=Sig;
for ii=1:iter
    
    Sig=extr_Sig1;
    Sigg=Sig;

    extr_Sig1(iii)=0;
    IF_est= ADTFD_IF_estimation_viterbi_modified(Sig, num,delta);
    for i=1:num
        
        IF=IF_est(i,:);
        
        Phase=2*pi*filter(1,[1 -1],IF);
        s_dechirp=exp(-1i*Phase);
        L=bw;
        %TF filtering for each sensor
        s1 = Sig.*(s_dechirp);
        s2=fftshift(fft(s1));
        s3=zeros(1,length(Sig));
        s3(WIN-L:WIN+L)=s2(WIN-L:WIN+L).*hamming(2*L+1)';
        s2(WIN-L:WIN+L)=0;
        extr_Sig=ifft(ifftshift(s3)).*conj(s_dechirp);
        s2=ifft(ifftshift(s2)).*conj(s_dechirp);
        Sig=s2;%-extr_Sig(iii);
        extr_Sig1(iii)=extr_Sig1(iii)+extr_Sig(iii);
              
        
    end
    if ii==2
        M1=mean(abs(extr_Sig1-Seg).^2)
        [tfd,orienttfd]=HTFD_new1(Sigg,2,50,84);
        
        figure; SetFigDef(16,9,'Times',20); tfsapl(real(extr_Sig1),tfd,'YLabel','Sample Number','Title','(c)', 'TFfontSize' , 20,'grayscale','on');
        %%%%%%%%%%%%%%%%
        t=0:127;
        figure
        SetFigDef(16,9,'Times',20);
        % set(gcf,'Color','w');
        plot(IF_est(:,1:4:end),t(1:4:end),'o',IF_O',t,'-','linewidth',2);
        ylabel('Sample Number');
        xlabel('Instantaneous Frequency (Hz)');
        title('(d)');
        axis([0 0.5 0 128]);
        
    end
    
    
end



for ii=1:2*iter
    Sig=extr_Sig1;
    extr_Sig1(iii)=0;
    
    for i=1:num
        
        IF=IF_est(i,:);
        
        Phase=2*pi*filter(1,[1 -1],IF);
        s_dechirp=exp(-1i*Phase);
        L=bw;
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
       
  
end
