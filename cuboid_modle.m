% close all
% clc
% plane1 = [0.1,0.2,0.3,0.4];
% s = rng;
% x = rand(1,1000);
% y = rand(1,1000);
% z = -(plane1(1)*x + plane1(2)*y + repmat(plane1(4), 1, length(x)));
% 
% figure
% scatter3(x,y,z,'.')
% view(-30,10)

%% simulate a cubiod
clc
clear
close all
%width,height and depth
width = 0.7;
height = 0.5;
depth = 0.3;

rotation = rotx(30)*roty(30)*roty(50);
translation =[0 1 1];

%x-y plane
plane_xy_z = repmat(depth/2,uint16(width*height*10000),1);
plane_xy_x = randi([-width * 500, width*500], uint16(width*height*10000),1)/1000;
plane_xy_y = randi([-height * 500, height*500], uint16(width*height*10000),1)/1000;
%x-z plane
plane_xz_y = repmat(height/2,uint16(width*depth*10000),1);
plane_xz_x = randi([-width * 500, width*500], uint16(width*depth*10000),1)/1000;
plane_xz_z = randi([-depth * 500, depth*500], uint16(width*depth*10000),1)/1000;
%y-z plane
plane_yz_x = repmat(width/2,uint16(depth*height*10000),1);
plane_yz_y = randi([-height * 500, height*500], uint16(depth*height*10000),1)/1000;
plane_yz_z = randi([-depth * 500, depth*500], uint16(depth*height*10000),1)/1000;

%draw cub
figure
hold on
% scatter3(plane_xy_x,plane_xy_y,plane_xy_z,'r.')
% scatter3(plane_xz_x,plane_xz_y,plane_xz_z,'g.')
% scatter3(plane_yz_x,plane_yz_y,plane_yz_z,'b.')


%cuboid cloud
cuboid(:,1) = [plane_xy_x;plane_xz_x;plane_yz_x];
cuboid(:,2) = [plane_xy_y;plane_xz_y;plane_yz_y];
cuboid(:,3) = [plane_xy_z;plane_xz_z;plane_yz_z];
%figure
scatter3(cuboid(:,1),cuboid(:,2),cuboid(:,3),'.')
view(-30,10)


cuboid_trans = cuboid*rotation +repmat(translation,size(cuboid,1),1);
% hold on
% scatter3(cuboid_trans(:,1),cuboid_trans(:,2),cuboid_trans(:,3),'r.')
center =[0,0,0];
z_axis = [0 0 1];
y_axis = [0 1 0];
x_axis = [1 0 0];
z_center = [center;z_axis];
y_center = [center;y_axis];
x_center = [center;x_axis];
plot3(x_center(:,1),x_center(:,2),x_center(:,3),'r');
plot3(y_center(:,1),y_center(:,2),y_center(:,3),'g');
plot3(z_center(:,1),z_center(:,2),z_center(:,3),'b');

%% wangyue algorithm
% find planes
points = cuboid_trans';
% points = loadpcd('../segment_pcd_result/cloud_cluster_15.pcd'); % box 2 9 11
[n1,b1,bestn1] = fitPlane(points,3,100);
p1 = points(:,bestn1);
% hold on;scatter3(points(1,bestn1),points(2,bestn1),points(3,bestn1),'r.');
points = points(:,~bestn1);
[n2,b2,bestn2] = fitPlane(points,3,100);
p2 = points(:,bestn2);
% hold on;scatter3(points(1,bestn2),points(2,bestn2),points(3,bestn2),'g.');
p3 = points(:,~bestn2);
pp{1} = p1;
pp{2} = p2;
pp{3} = p3;
n{1} =n1;
n{2} =n2;
n{3} = cross(n1,n2);
axis equal
plane_size = [length(p1), length(p2), length(p3)];
[plane_size_sort, plane_size_index] = sort(plane_size, 'descend');
p1 = pp{plane_size_index(1)};
p2 = pp{plane_size_index(2)};
p3 = pp{plane_size_index(3)};
n1 = n{plane_size_index(1)};
n2 = n{plane_size_index(2)};
n3 = n{plane_size_index(3)};
%% 
% p1 = pp{1};
% p2 = pp{3};
% p3 = pp{2};
% n1 = n{1};
% n2 = n{3};
% n3 = n{2};
% standard align
% n1 aligned with z
t = acos(n1(3));
u = cross(n1,[0 0 1].');
u = u./norm(u);
R = cos(t)*eye(3)+sin(t)*[0 -u(3) u(2);u(3) 0 -u(1);-u(2) u(1) 0]+(1-cos(t))*[u(1)*u(1) u(1)*u(2) u(1)*u(3);u(2)*u(1) u(2)*u(2) u(2)*u(3);u(3)*u(1) u(3)*u(2) u(3)*u(3)];
% hold on;scatter3(p1(1,:),p1(2,:),p1(3,:),'r.');
tp1 = R*p1;
%T = -[0;0;mean(tp1(3,:))];
T = -[mean(tp1(1,:));mean(tp1(2,:));0];
tp1 = R*p1+repmat(T,1,size(p1,2));
tp2 = R*p2+repmat(T,1,size(p2,2));
% hold on;scatter3(tp1(1,:),tp1(2,:),tp1(3,:),'r.');
% hold on;scatter3(tp2(1,:),tp2(2,:),tp2(3,:),'g.');

% projected n2 aligned with y
n2n = R*n2;
t = pi/2 + atan2(n2n(2),n2n(1));
c = cos(t);s = sin(t);
R2 = [c -s 0;s c 0;0 0 1];
%R2 = R2';
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

T2 = -[mean(tp2(1,:));0;mean(tp2(3,:))];
% tf = T2;
tp1 = tp1 + repmat(T2,1,size(tp1,2));
tp2 = tp2 + repmat(T2, 1, size(tp2,2));
tf = tf+T2;
p = [tp1 tp2];
% scatter3(tp1(1,:),tp1(2,:),tp1(3,:),'r*');
% scatter3(tp2(1,:),tp2(2,:),tp2(3,:),'g*');
w = quantile(p(1,:),0.999)-quantile(p(1,:),0.001);
l = quantile(p(2,:),0.999)-quantile(p(2,:),0.001);
h = quantile(p(3,:),0.999)-quantile(p(3,:),0.001);
%h = quantile(abs([p(3,:) p(3,:)]),0.999);
%op = [min(p(1,:));min(p(2,:));min(p(3,:))];
op = [-w;-height;-depth];
op = 0.5*op;
box = [op op+[w 0 0]' op+[0 l 0]' op+[w l 0]' op+[0 0 h]' op+[w 0 h]' op+[0 l h]' op+[w l h]'];
tf = -Rf'*tf;
box = Rf'*box + repmat(tf,1,8);
% draw
drawbox(box,'k')

%% my algorithm
% caculate the maxmize plane width and height
rotation1=[n3(1),n3(2),n3(3);n2(1),n2(2),n2(3);n1(1),n1(2),n1(3)];
M = zeros(4,3);
D = zeros(4,1);

center_major = [mean(p1(1,:));mean(p1(2,:));mean(p1(3,:))];
M(1,1) = n1(2);M(1,2) = -n1(1);
M(2,2) = n1(3);M(2,3) = -n1(2);
D(1,1) = M(1,:)*center_major;
D(2,1) = M(2,:)*center_major;

center_middle = [mean(p2(1,:));mean(p2(2,:));mean(p2(3,:))];
M(3,1) = n2(2);M(3,2) = -n2(1);
M(4,2) = n2(3);M(4,3) = -n2(2);
D(3,1) = M(3,:)*center_middle;
D(4,1) = M(4,:)*center_middle;


w = width;
l = height;
h = depth;
op = [-w;-height;-depth];
op = 0.5*op;
box = [op op+[w 0 0]' op+[0 l 0]' op+[w l 0]' op+[0 0 h]' op+[w 0 h]' op+[0 l h]' op+[w l h]'];

translation1 = inv((M'*M))*M'*D;
% 
 box = rotation1'*box + repmat(translation1,1,8);
% draw
 figure
hold on

scatter3(p1(1,:),p1(2,:),p1(3,:),'r.');
scatter3(p2(1,:),p2(2,:),p2(3,:),'g.');
scatter3(p3(1,:),p3(2,:),p3(3,:),'b.');
drawbox(box,'b')
center =[0,0,0];
z_axis = [0 0 1];
y_axis = [0 1 0];
x_axis = [1 0 0];
z_center = [center;z_axis];
y_center = [center;y_axis];
x_center = [center;x_axis];
plot3(x_center(:,1),x_center(:,2),x_center(:,3),'r');
plot3(y_center(:,1),y_center(:,2),y_center(:,3),'g');
plot3(z_center(:,1),z_center(:,2),z_center(:,3),'b');
scatter3(cuboid(:,1),cuboid(:,2),cuboid(:,3),'.')
% scatter3(plane_xy_x,plane_xy_y,plane_xy_z,'r.')
% scatter3(plane_xz_x,plane_xz_y,plane_xz_z,'g.')
% scatter3(plane_yz_x,plane_yz_y,plane_yz_z,'b.')
view(-30,10)
%% error
rotation_error1 = norm((abs(rotation) - abs(Rf)).^2)
rotation_error2 = norm((abs(rotation) - abs(rotation1)).^2)

trans_error1 = norm((translation' - tf).^2)
trans_error2 = norm((translation' - translation1).^2)



















