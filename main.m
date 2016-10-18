% see 2 planes
close all;clc;clear
points = loadpcd('./segment_pcd_result/cloud_cluster_2.pcd'); % box 2 9 11
fullpoints = points;

% find planes
[n1,b1,bestn1] = fitPlane(points,3,100);
p1 = points(:,bestn1);
% hold on;scatter3(points(1,bestn1),points(2,bestn1),points(3,bestn1),'r.');
points = points(:,~bestn1);
[n2,b2,bestn2] = fitPlane(points,3,100);
p2 = points(:,bestn2);
% hold on;scatter3(points(1,bestn2),points(2,bestn2),points(3,bestn2),'g.');
p3 = points(:,~bestn2);
axis equal

% standard align
% n1 aligned with z
t = acos(n1(3));
u = cross(n1,[0 0 1].');
u = u./norm(u);
R = cos(t)*eye(3)+sin(t)*[0 -u(3) u(2);u(3) 0 -u(1);-u(2) u(1) 0]+(1-cos(t))*[u(1)*u(1) u(1)*u(2) u(1)*u(3);u(2)*u(1) u(2)*u(2) u(2)*u(3);u(3)*u(1) u(3)*u(2) u(3)*u(3)];
% hold on;scatter3(p1(1,:),p1(2,:),p1(3,:),'r.');
tp1 = R*p1;
T = -[0;0;mean(tp1(3,:))];
tp1 = R*p1+repmat(T,1,size(p1,2));
tp2 = R*p2+repmat(T,1,size(p2,2));
% hold on;scatter3(tp1(1,:),tp1(2,:),tp1(3,:),'r.');
% hold on;scatter3(tp2(1,:),tp2(2,:),tp2(3,:),'g.');

% projected n2 aligned with y
n2n = R*n2;
t = atan2(n2n(2),n2n(1));
c = cos(t);s = sin(t);
R2 = [c -s 0;s c 0;0 0 1];
tp1 = R2'*tp1;
tp2 = R2'*tp2;
% hold on;scatter3(tp1(1,:),tp1(2,:),tp1(3,:),'r.');
% hold on;scatter3(tp2(1,:),tp2(2,:),tp2(3,:),'g.');

% total projection
Rf = R2'*R;
tf = R2'*T;
% fullpoints = Rf*fullpoints+repmat(tf,1,size(fullpoints,2));
hold on;
scatter3(p1(1,:),p1(2,:),p1(3,:),'r.');
scatter3(p2(1,:),p2(2,:),p2(3,:),'g.');
scatter3(p3(1,:),p3(2,:),p3(3,:),'b.');
%scatter3(fullpoints(1,:),fullpoints(2,:),fullpoints(3,:),'g.');

% get recangle
p = [tp1 tp2];
w = quantile(p(1,:),0.999)-quantile(p(1,:),0.001);
l = quantile(p(2,:),0.999)-quantile(p(2,:),0.001);
h = quantile(p(3,:),0.999)-quantile(p(3,:),0.001);
%h = quantile(abs([p(3,:) p(3,:)]),0.999);
op = [min(p(1,:));min(p(2,:));min(p(3,:))];
box = [op op+[w 0 0]' op+[0 l 0]' op+[w l 0]' op+[0 0 h]' op+[w 0 h]' op+[0 l h]' op+[w l h]'];
box = Rf'*(box-repmat(tf,1,8));
% draw
drawbox(box,'k')

% refinement



