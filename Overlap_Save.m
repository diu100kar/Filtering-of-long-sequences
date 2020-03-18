function answ=ovsave(x,h,L)
    lx=length(x);
    M=length(h);
    N=M+L-1;
    H=fft(h,N);
    x2=zeros(1,ceil(length(x)/L)*L);
    x2(1:length(x))=x;
    x2=[zeros(1,M-1),x2,zeros(1,M-1)];
    answ=[];
    for i=1:L:lx+M-1
        tempfft=fft(x2(1,i:i+N-1),N);
        y=ifft(H.*tempfft,N);
        answ=[answ,y(1,M:end)];
    end
    answ=round(answ(1,1:lx+M-1));
end
