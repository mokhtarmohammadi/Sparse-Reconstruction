function [ S_ADTFD ] = sparse_adtfd( adtfd,num )

M=length(adtfd);
n=0:M-1;

for ii=0:M-1
    A(ii+1,:)=exp(1i*ii*2*pi*n/M)/M;
end

%adtfd11=ifft(HTFD_new1(x1,2,30,64),[],1);

adtfd=ifft(adtfd,[],1);
LL=num;
%I=quadtfd(x,length(x)-1,1,'specx',85,'hamm',128);
T=0.1;
S_ADTFD=zeros(size(adtfd));
for nn=1:M
                aa=adtfd(:,nn);
              
        index=find(abs(aa)>=T*max(abs(aa)));
        A1=A(index,:);
        if length(index)<LL
            LL=length(index);
        end
            S_ADTFD(:,nn)=real((OMP( A1, adtfd(index,nn), LL )));
    

end

