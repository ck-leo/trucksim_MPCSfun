function [angle] = NormalizeAngle(angle)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
a = mod(angle+pi,2*pi);
if(a<0)
    a = a+2*pi;
end
angle = a-pi;
end

