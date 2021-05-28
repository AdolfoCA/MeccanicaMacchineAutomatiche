clc
clear
close all
 
global l d
l=1;
d=0.5;

Q=[2*l;1*l;pi/2+0.01];

P=kin_dir_pos(Q)


plot([Q(1) Q(1)+l*cos(Q(3))],[0 l*sin(Q(3))])
hold on
plot([Q(1)+l*cos(Q(3)) Q(1)+l*cos(Q(3))+l*cos(P(3)+pi/2)],[l*sin(Q(3)) l*sin(Q(3))+l*sin(P(3)+pi/2)])
plot([Q(1)+l*cos(Q(3))+l*cos(P(3)+pi/2) 0],[l*sin(Q(3))+l*sin(P(3)+pi/2) Q(2)])

plot([Q(1)+l*cos(Q(3))+l/2*cos(P(3)+pi/2) P(1)],[l*sin(Q(3))+l/2*sin(P(3)+pi/2) P(2)])
plot(P(1),P(2),'o')
xlim([-0.5 2.5])
ylim([-0.5 2.5])



