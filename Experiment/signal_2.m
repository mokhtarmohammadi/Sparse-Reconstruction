function [x y iii iiii]=signal_2(SampFreq,M)
                        
% SampFreq = 128; sample
fs=SampFreq;
t = 0:1/SampFreq:1-1/SampFreq;
N=SampFreq;
n = 0:127;
% M=64;
%M=32;
iter=10;
Sig1 = exp(1i*(2*pi*(0.1*n.^3)/(128*128))+1i*(2*pi*(0.05*n))+1i*(2*pi*(0.05*n.^2/128))); %300t�???150t
Sig2= exp(1i*(2*pi*(0.1*n.^3)/(128*128))+1i*(2*pi*(0.15*n))+1i*(2*pi*(0.05*n.^2/128)));
%  Sig2= exp(1i*(-2*pi*(0.1*n.^3)/(128*128))+1i*(2*pi*(0.45*n))-1i*(2*pi*(0.05*n.^2/128)));
IF_O(:,1)=(.3/(128*128))*n.^2+.05+(.1/128)*n;
IF_O(:,2)=(.3/(128*128))*n.^2+.15+(.1/128)*n;
Sig=Sig1+Sig2;
seg=Sig;
p = randperm(N);
switch M
    case 16
p(1:M)=[10:14 60:64 90:95];
    case 32
        p(1:M)=[10:19 60:69 90:101];
    case 48
        p(1:M)=[10:24 60:74  90:107];%84:84+7 100:107];
    case 64
p(1:M)=[9:32 60:83  90:105 ];
    otherwise
end
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