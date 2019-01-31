function [extr_Sig1,fidexmult, tfdv] = extridge_mult_tffiltering_tfds(Sig,TFD, SampFreq, num, delta,bw,iiii,iii)
fs=SampFreq;
if (isreal(Sig))
Sig = hilbert(Sig);
end

orderamp = round(bw*length(Sig)/SampFreq);%Fourier order for characterizing signal amplitude
fidexmult = zeros(num,length(Sig));
tfdv = zeros(num,length(Sig));
extr_Sig1=zeros(1,length(Sig));
extr_Sig1=Sig;
%               am1 =wvd1(Sig);
% %% computation of theta
% theta = -90:2:90;
% [R,xp] = radon(abs(am1),theta);
% pos=find(xp==0);
% % R(pos,:)=R(pos,:)/sum(R(pos,:));
% Rad_0=Dtrend(R(pos,:),10);
% R_trend = R(pos,:)-Rad_0';
% Rad_0 = Rad_0/max(Rad_0);
% % Rad_0=R(pos,:)/max(R(pos,:));
% [ peak,loc ] = Detection_pic(Rad_0);
% pos_pic = find(peak>0.3);
% angles = 180-theta(loc(pos_pic));

for i = 1:num
%%%%%%%%%%%%%%%%%%

f = linspace(0,SampFreq/2,length(extr_Sig1));
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
 Orient=orienttfd;
Spec(and(80<Orient,Orient<100))=0;
I2=Spec;
         case 'F_ADTFD'
             
%                 [Spec2,orienttfd2]=HTFD_neww(Sig,2,15,84);
%         [Spec1,orienttfd1]=HTFD_neww(Sig,2,30,84);
          [Spec2,orienttfd2]=F_ADTFD(Sig,2,15,84,angles);
          [Spec1,orienttfd1]=F_ADTFD(Sig,2,30,84,angles);
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
 Orient=orienttfd;
Spec(and(80<Orient,Orient<100))=0;
I2=Spec;
end
I2(I2<0)=0;
Spec=I2;

switch TFD
case 'ADTFD'
c = findridges_new1(Spec,orienttfd,2*delta);
 case 'F_ADTFD'
           c = findridges_new1(Spec,orienttfd,2*delta);
otherwise

           c = findridges(Spec,2*delta);
end                 
IF=(length(Sig)-c)/(2*length(Sig));
                 IF=(c)/(2*length(Sig));
                Phase=2*pi*filter(1,[1 -1],IF);
                s_dechirp=exp(-1i*Phase);
                L=4;
                %TF filtering for each sensor
                s1 = Sig.*(s_dechirp);
                s2=fftshift(fft(s1));
                s3=zeros(1,length(Sig));
                s3(128-L:128+L)=s2(128-L:128+L);
                s2(128-L:128+L)=0;
                extr_Sig=ifft(ifftshift(s3)).*conj(s_dechirp);
                s2=ifft(ifftshift(s2)).*conj(s_dechirp);
                Sig=s2;%-extr_Sig(iii);
               extr_Sig1(iii)=extr_Sig1(iii)+extr_Sig(iii);
 % remove the extracted signal component so that other ridge curves with smaller energies can be extracted in the subsequent iterations
end
% extr_Sigtemp=extr_Sig1*mean(abs(Sig(iiii)))/mean(abs(extr_Sig1(iiii)));
% %extr_Sig1=extr_Sig1*exp(1i*angle(sum(extr_Sig1.*conj(Sig))));
% extr_Sig=extr_Sigtemp*exp(1i*angle(sum(extr_Sigtemp(iiii).*conj(Sig(iiii)))));
% save extridge Spec1 Spec2 Spec3 Sig1 Sig2 Sig3;