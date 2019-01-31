function IF = findridges_new_viterbi_adtfd(Spec,Specangle)
Path=zeros(size(Spec));
Pathweight=zeros(1,length(Spec));
Specangle(Specangle>90)=-180+Specangle(Specangle>90);
Path(:,1)=1:length(Path);
for i=2:length(Spec)
    for ii=1:length(Spec) % all members of
        Pathw=zeros(1,length(Spec));
        % For transition from ii to j.
        
        for jj=1:length(Spec)
            rr=abs(jj-ii)-2;
            if rr<0
                rr=0;
                
            end
            
%             if angg<5
%                 angg=0;
%             end
                        angg=abs(Specangle(jj,i-1)-Specangle(ii,i));
            if angg<5
                angg=0;
            else
                angg=(angg-5)*30;
                      %          angg=angg*5;

            end
%             if i>2
%             angg=Path(jj,i-1)-jj-(jj-ii);
%             if (Path(jj,i-1)-jj)*(jj-ii)>0
%                 angg=50*angg;
%             else
%                 angg=0;
%             end
%             else
%                 angg=0;
%             end
            Pathw(jj)=rank(Spec(:,i),ii)+Pathweight(jj)+1000*rr+angg;% weight of ii +
                   %     Pathw(jj)=rank(Spec(:,i),ii)+Pathweight(jj)+5*rr+angg;% weight of ii +

        end
        [value,index]=min(Pathw);
        Pathweight1(ii)=value;
        Path(ii,i)=index;
    end
    Pathweight=Pathweight1;
    
end
IF(1,length(Spec))=0;
[~,index]=min(Pathweight);
for i=length(Spec):-1:1
    IF(i)=index;
    index=Path(index,i);
end
% Code for decoding
end
function a= rank(A,b)
As=sort(A,'descend');
ind=find(As==A(b));
%if A(b)~=0
a=ind(1)-1;
%else
 %   a=128;

%end
end