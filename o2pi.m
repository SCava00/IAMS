function [th]=o2pi(th)
while(th>=2*pi)
    th=th-2*pi;
end
while(th<0)
    th=th+2*pi;
end