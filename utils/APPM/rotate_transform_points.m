% generate three random points
Px1 = rand(1)*10;
Py1 = rand(1)*10;

Px2 = rand(1)*10;
Py2 = rand(1)*10;

Px3 = rand(1)*10;
Py3 = rand(1)*10;

% random point matrix
P1 = [Px1; Py1];
P2 = [Px2; Py2];
P3 = [Px3; Py3];

% generate random rotation matrix about x-axis

r = 90;  % rotation value

R = [cos(r) -sin(r); sin(r) cos(r)];

% rotated points
M1 = R*P1;
M2 = R*P2;
M3 = R*P3;

% translation matrix
T = [0.5; 2]

% plot points

figure(1)
hold on

% random plotted points
plot([Px1 Px2 Px3],[Py1 Py2 Py3],'r')

% rotated+transformed plotted points
plot([M1(1)+T(1) M2(1)+T(1) M3(1)+T(1)],[M1(2)+T(2) M2(2)+T(2) M3(2)+T(2)],'b')


hold off