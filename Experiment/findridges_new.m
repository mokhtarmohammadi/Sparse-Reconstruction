function index = findridges_new(Spec,Specangle,delta)
%Ridge detection algorithm,i.e., Algorithm 1 in paper:Separation of Overlapped Non-Stationary Signals by Ridge Path Regrouping and Intrinsic Chirp Component Decomposition
% ,IEEE Sensors journal,2017.
% The algorithm is originally introduced in the paper: Algorithms for blind components separation and extraction from the time-frequency distribution of their mixture.
% EURASIP Journal on Advances in Signal Processing, 2004.


%Spec£ºTime-Frequency distribution of the signal
%delta£ºmaximum allowable frequency variation between two consecutive points
%index£ºThe obtained frequency indexs at each time instant
Spec = abs(Spec);
[M,N] = size(Spec);
index = zeros(1,N);
index1=index;

[fmax,tmax] = find(Spec == max(Spec(:)));
fmax = fmax(1);
tmax = tmax(1);
index(tmax) = fmax;
%Specangle(Specangle>90)=-180+Specangle(Specangle>90);
%Specangle=Specangle-90;
Specangle(fmax,tmax)=Specangle(fmax,tmax);
theta=Specangle(fmax,tmax);
%theta=180-theta;
T=15;
f0 = fmax;
jj=-1;
for j = (min(tmax+1,N)):N
    f0=round(f0);
    
    low = max(1,f0-delta);
    up = min(M,f0+delta);
    [~,f0] = max(Spec(low:up,j));
    %    [~,f0] = min(abs(Specangle(low:up,j)-theta));
    
    f0 = f0 + low - 1;
    
    if abs(Specangle(round(f0),j)-theta)>T
    %if min(abs(Specangle(round(f0),j)-theta),abs(180-Specangle(round(f0),j)-theta))>T
        if jj==-1
            jj=j-1;
        end
%         f0=index(j-1)-tan(pi*theta/180);
        f0=index(jj)-(j-jj)*tan(pi*theta/180);
          %      f0=index(jj)-(j-jj)/tan(pi*theta/180);
%           fo=index(j)-index(max(j-8,1))/10;
    else
        jj=-1;
        %if Specangle(round(f0),j)~=-10
            theta=0.6*theta+0.4*Specangle(round(f0),j);
            %        theta=Specangle(round(f0),j);
            
       % end
        
    end
    %  theta=0.9*theta+0.1*Specangle(round(f0),j);
    
    
    if f0<1
        f0=1;
    end
    if f0>255
        f0=255;
    end
    index(j) = f0;
%     index1(j)=theta;
    %f0
end
theta=Specangle(fmax,tmax);
jj=-1;

f1 = fmax;
for j = (max(1,tmax-1)):-1:1
    f1=round(f1);
    low = max(1,f1-delta);
    up = min(M,f1+delta);
    [~,f1] = max(Spec(low:up,j));
    f1 = f1 + low - 1;
    %if abs(Specangle(round(f1),j)-Specangle(round(index(j+1)),j+1))>T
    if abs(Specangle(round(f1),j)-theta)>T
   % if min(abs(Specangle(round(f1),j)-theta),abs(180-Specangle(round(f1),j)-theta))>T   
        %  f1=index(j+1)+tan(pi*theta/180);
        if jj==-1
            jj=j+1;
        end
        
        f1=index(jj)-(j-jj)*tan(pi*theta/180);
%                        f1=index(jj)-(j-jj)/tan(pi*theta/180);
% f1=index(j)-index(max(j-8,1))/10;
    else
        jj=-1;
       % if Specangle(round(f1),j)~=-10
            
            theta=0.6*theta+0.4*Specangle(round(f1),j);
            %    theta=Specangle(round(f1),j);
            
       % end
    end
    %   theta=0.9*theta+0.1*Specangle(round(f0),j);
    
    
    if f1>255
        f1=255;
    end
    if f1<1
        f1=1;
    end
    index(j) = f1;
end
index=round(index);
end

