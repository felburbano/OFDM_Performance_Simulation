function ecualized = ZF_Equalize(N,r,u)
    o=fft(r,N);    
    ecualized=u./o;
end