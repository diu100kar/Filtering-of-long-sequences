function answ=ovadd(x,h,L)
    lx=length(x);
    M=length(h);
    N=M+L-1;
    h=[h,zeros(1,L-1)];
    H=fft(h,N);
    x2=zeros(1,ceil(length(x)/L)*L);
    x2(1:length(x))=x;
    answ=zeros(1,length(x2)+M-1);
    for i=1:L:length(x)
        xtemp=[x2(1,i:i+L-1),zeros(1,M-1)];
        tempfft=fft(xtemp,N);
        tempY=tempfft.*H;
        y=ifft(tempY,N);
        answ(1,i:i+N-1)=answ(1,i:i+N-1)+y;
    end
    answ=round(answ(1,1:M+lx-1));
end
