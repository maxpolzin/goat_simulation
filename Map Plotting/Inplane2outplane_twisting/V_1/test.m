for ii=1:100
    for jj=1:10
        F1=figure(1)
        x=linspace(0,ii*pi,ii*120)
        y1=sin(x)
        plot(x,y1)
        pause(0.5)
    end
end