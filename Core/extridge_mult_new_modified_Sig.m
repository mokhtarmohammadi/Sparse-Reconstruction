function [fidexmult Spec]  = extridge_mult_new_modified_Sig(Sig, num, delta)
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

for i = 1:num
%Spec=tfr_stft_high(Sig);

                [Spec2,orienttfd2]=HTFD_neww(Sig,2,15,84);
        [Spec1,orienttfd1]=HTFD_neww(Sig,2,30,84);

Spec=min(Spec1,Spec2);
for ii=1:length(Spec1)
    for jj=1:length(Spec2)
            value=min(Spec1(ii,jj),Spec2(ii,jj));
            Spec(ii,jj)=value;
            if Spec1(ii,jj)==value
            orienttfd(ii,jj)=orienttfd1(ii,jj);
            else
            orienttfd(ii,jj)=orienttfd2(ii,jj);

            end
            
         
    end
end
% 
    %    [Spec,orienttfd1]=HTFD_neww(Sig,2,30,84);

c = findridges_new1(Spec,orienttfd,2*delta);
%c = findridges_new(Spec,orienttfd,2*delta);

%sc=findridges_new(Spec,orienttfd,delta);
%c=findridges(Spec,2*delta);

 IF=(c)/(2*length(Sig));

                Phase=2*pi*filter(1,[1 -1],IF);
                s_dechirp=exp(-1i*Phase);

                
                %im_label2=bwmorph(im_label2,'dilate',3);
                
                % For each sensor do the following steps
                
                L=2;
                %TF filtering for each sensor
                s1 = Sig.*(s_dechirp);
                s2=fftshift(fft(s1));
                s3=zeros(1,length(Sig));
                s3(128-L:128+L)=s2(128-L:128+L);
                s2(128-L:128+L)=0;
                extr_Sig=ifft(ifftshift(s3)).*conj(s_dechirp);
                s2=ifft(ifftshift(s2)).*conj(s_dechirp);
                
                %Sig(iii)=Sig(iii)-extr_Sig(iii);
                Sig=s2;%-extr_Sig(iii);

           % extr_Sig1=extr_Sig1+extr_Sig;
         
%  figure; imagesc(Spec);hold on; quiver(1:256,1:256,cos(orienttfd*pi/180),sin(orienttfd*pi/180));
fidexmult(i,:) = c;
extr_Sigtemp=extr_Sig*mean(abs(Sig))/mean(abs(extr_Sig));
%extr_Sig1=extr_Sig1*exp(1i*angle(sum(extr_Sig1.*conj(Sig))));
extr_Sig=extr_Sigtemp*exp(1i*angle(sum(extr_Sigtemp.*conj(Sig))));
SpecN=HTFD_neww(extr_Sig,2,30,84);
SpecT(i,:,:)=SpecN/max(SpecN(:));

end
Spec=SpecT(1,:,:)+SpecT(2,:,:)+SpecT(3,:,:)+SpecT(4,:,:);
 
end