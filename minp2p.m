function [th]=minp2p(th)
while(th<-pi)
    th=th+2*pi;
end
while(th>pi)
    th=th-2*pi;
end