function [y,it,An,MSE]=GradRec(y,Nxx,PR,Delta,Nit,x0)
% Gradient based reconstruction of missing samples
% Milos Dakovic, 29.3.2014.
%
% Usage:
%       y = GradRec(...)
%       [y,it] = GradRec(...)
%       [y,it,An] = GradRec(...)
%       [y,it,An,MSE] = GradRec(...)
%       ... = GradRec(x,Nxx)
%       ... = GradRec(x,Nxx,PR)
%       ... = GradRec(x,Nxx,{@transform,PR})
%       ... = GradRec(x,Nxx,PR,Delta)
%       ... = GradRec(x,Nxx,PR,Delta,Nit)
%       ... = GradRec(x,Nxx,PR,Delta,Nit,x0)
%
% Inputs:
%       x         - signal - row vector
%       Nxx       - positions of missing samples
%       PR        - target precision (in dB) or cell array with precision and transform function
%                   default precision is 100dB and default transform is fft
%                   examples: 150, {@dct}, {@fft,150}, {@UserDefinedTransform,80}
%       Delta     - gradient algorithm parameter
%                   default Delta is maximal amplitude of available samples
%       Nit       - maximal number of iterations per cycle
%                   default value is signal length
%       x0        - original signal, used for MSE calculations only
%                   if signal is not provided the MSE is not calculated
%
% Outputs:
%       y       - reconstructed signal
%       it      - number of iterations
%       An      - Angles between sucessive gradients
%       MSE     - MSE for each iteration
%

N=length(y);
K=length(Nxx);

% Input arguments
if nargin<3
    PR=100;
    Trans=@fft;
else
    if iscell(PR)
        Trans=PR{1};
        if length(PR)>1
            PR=PR{2};
        else
            PR=100;
        end
    else
        Trans=@fft;
    end
end
if nargin<4
    yy=y;yy(Nxx)=[];
    Delta=max(abs(yy));
end
if nargin<5
    Nit=N;
end


MSE=[];
An=[];
it=0;

% Form matrix D
D=zeros(N,K);
for k=1:K
    D(Nxx(k),k)=1;
end
D=feval(Trans,D)';

% Gradient algorithm
Gr=zeros(1,N);
y0=y;
for p=1:35
    Ind=3;
    DD=D*Delta;
    CC=2*Delta/(2*Delta*N);
    
    for k=1:Nit
        it=it+1;
        Y=repmat(feval(Trans,y).',1,K);
        
        gg=(sum(abs(Y+DD'))-sum(abs(Y-DD')));
        Gr(Nxx)=gg;    
        y=y-CC*Gr;
        
        % Calculate MSE if original signal is provided
        if nargin>5 
            MSE(end+1)=mean(abs(x0(Nxx)-y(Nxx)).^2);
        end
        
        % Stoping criterion
        Gri=sqrt(sum(Gr.^2));
        if Gri==0 
            % Stop if gradijent is null vector
            Ind=1;
            break
        end
        if k>1
            AAn=acos(sum(Gp.*Gr)/(Gpi*Gri))*180/pi;
            if nargout>2
                An(it)=AAn;
            end
            if AAn>170
                % Stop if angle between two successive gradients is higher that 170 degrees
                Ind=0;
                break
            end
        else
            if nargout>2
                % Angle is undefined in the first iteration
                An(it)=NaN;
            end           
        end        
        Gp=Gr; Gpi=Gri;
    end
    
    if Ind==0
        Delta=Delta/sqrt(10);
        TT=10*(log10(mean(abs(y(Nxx)).^2))-log10(mean(abs(y0(Nxx)-y(Nxx)).^2)));
        if TT>PR
            % Precision is achieved
            break
        end
        y0=y;
    else
%        warning('MATLAB:GradRec','The algorithm is stopped before required precision is reached.')
        break
    end
end


