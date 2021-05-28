clc
clear
close all
 
global l d
l=1;
d=0.5;

Q=[1*l;1.5*l;pi/4];

P=kin_dir_pos(Q)

plot([0 Q(1)],[0 0],'r')
hold on
plot([Q(1) Q(1)+l*cos(Q(3))],[0 l*sin(Q(3))])
plot([Q(1)+l*cos(Q(3)) Q(1)+l*cos(Q(3))+l*cos(P(3)+pi/2)],[l*sin(Q(3)) l*sin(Q(3))+l*sin(P(3)+pi/2)])