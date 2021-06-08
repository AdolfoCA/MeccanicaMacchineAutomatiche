close all

% qui ci vanno i parametri reali
l=0.99; d=0.49;

% Utilizzo i dati salvati dalla simulazione
Q=out.simout.Data(:,:);
P=out.simout1.Data(:,:);

% Grafica
for i=1:10:length(Q)
figure(1)
plot([0 0],[0 2],'-sk','LineWidth',2)   %pattino
hold on
plot([0 2],[0 0],'-sk','LineWidth',2)   %pattino
plot([Q(i,1) Q(i,1)+l*cos(Q(i,3))],[0 l*sin(Q(i,3))],'-o')
plot([Q(i,1)+l*cos(Q(i,3)) Q(i,1)+l*cos(Q(i,3))+l*cos(P(i,3)+pi/2)],[l*sin(Q(i,3)) l*sin(Q(i,3))+l*sin(P(i,3)+pi/2)],'-o')
plot([Q(i,1)+l*cos(Q(i,3))+l*cos(P(i,3)+pi/2) 0],[l*sin(Q(i,3))+l*sin(P(i,3)+pi/2) Q(i,2)],'-o')
plot([Q(i,1)+l*cos(Q(i,3))+l/2*cos(P(i,3)+pi/2) P(i,1)],[l*sin(Q(i,3))+l/2*sin(P(i,3)+pi/2) P(i,2)],'-')

plot(P(i,1),P(i,2),'.','MarkerSize',20)
plot(Q(i,1),0,'sb','MarkerSize',40)
%il punto rosso indica la presenza di un motore sul giunto
plot(Q(i,1),0,'.','MarkerSize',40) 
plot(0,Q(i,2),'sb','MarkerSize',40)
xlim([-0.1 2])
ylim([-0.1 2])
axis equal
grid on
drawnow
hold off
end
