%Matlab reading csv and 3d animation plotting
clc;clear;close;set(0,'DefaultFigureWindowStyle','docked');
duration = 10;

positions = readmatrix('3.5BodySimResults.csv');
x = 1:3:10;
y = 2:3:11;
z = 3:3:12;

X = positions(:,x);
Y = positions(:,y);
Z = positions(:,z);

psl = 3e8;
lims = [-psl psl];
t = {'Sun','Earth','Mars','Ship'};
c = [[0.9290 0.4587 0.1250];
     [0 0.4470 0.7410];
     [0.6350 0.0780 0.1840];
     [0.4660 0.6740 0.1880]];
ms = [15 1 1 1];

% Unanimated plot of the orbits
figure(1); hold on;
for i=1:length(x(1,:))
    plot3(X(:,i),Y(:,i),Z(:,i),'.','Color',c(i,:), ...
          'MarkerSize',ms(i))
end
axis([lims lims lims]); xlabel('X'); ylabel('Y'); zlabel('Z');
grid on; grid minor; view(-45,45);

fig2 = figure(2);
axis([lims lims lims]); xlabel('X'); ylabel('Y'); zlabel('Z');
grid on; grid minor; hold on; 
n = length(positions(:,1));

for i=1:n
    for j=1:length(x(1,:))
        plot3(X(i,j),Y(i,j),Z(i,j),'.','Color',c(j,:), ...
              'MarkerSize',ms(j),'LineWidth',1)
    end
%   text(X(i,:),Y(i,:),Z(i,:),t)
    view(-45+0.1*i,45);

    frame(i) = getframe(fig2);
end
V = VideoWriter('3.5BodyOrbitsSpin.mp4','MPEG-4');
V.FrameRate = n/duration;
V.Quality = 50;
open(V);
writeVideo(V,frame);
close(V);