function [angle] = NormalizeAngle(angle)
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
a = mod(angle+pi,2*pi);
if(a<0)
    a = a+2*pi;
end
angle = a-pi;
end

