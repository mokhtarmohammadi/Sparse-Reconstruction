function SM=MySTFT(x)
L=5;
[M N]=size(x);
if M==1
    win=[0 sqrt(hanning(length(x)-1))']; Ew=sum(win.^2)/length(win);
    STFT=fftshift(fft(fftshift(x.*win))).'; %STFT
    stft=STFT/sqrt(Ew)/length(win); %Normalization for the plot
    SM=SM_calc(stft,L);
else
         for n=[1:N]
        X=x(:,n)';
        win=[0 sqrt(hanning(length(x)-1))']; Ew=sum(win.^2)/length(win);
        STFT=fftshift(fft(fftshift(X.*win))).'; %STFT
        stft=STFT/sqrt(Ew)/length(win); %Normalization for the plot
        SM(:,n)=SM_calc(stft,L);
            end
end
SM=SM';