clear all;
close all;

%last updated: 8.4.2020

%given information
d = [.340 0 .400 0 .400 0 .126]';
alpha = [-pi/2 pi/2 pi/2 -pi/2 -pi/2 pi/2 0]';
a = zeros(7,1);

%frequency = 100hz
%Gait Cycle = 140 frames = 1.4s
t = 1:1:140;

%approximate link diameter:
width = .068*2; %m

AverageCurvemat = load('AverageCurve.mat');
AverageCurve = AverageCurvemat.AverageCurve/1000; %m
RAmat = load('RA.mat');
RA = RAmat.RA;
% EulerRAmat = load('EulerRA.mat');
% EulerRA = EulerRAmat.EulerRA;
%initial angular configuration
qq = zeros(7,1);
q = zeros(length(qq), length(t));
q(:,1) = qq;

%Kuka limitations lim(joint) = (+-Range of motion (radians), max velocity (radians/s))
lim = zeros(7,2);
lim(1,1:2) = [170, 98]*pi/180;  
lim(2,1:2) = [120, 98]*pi/180;  
lim(3,1:2) = [170, 100]*pi/180;  
lim(4,1:2) = [120, 130]*pi/180;  
lim(5,1:2) = [170, 140]*pi/180;  
lim(6,1:2) = [120, 180]*pi/180;  
lim(7,1:2) = [170, 180]*pi/180;  

%Calculating Transformation Matricies
TB0 = [0 0 1 0; 1 0 0 0; 0 1 0 0; 0 0 0 1];
T07 = cell(140,1);
TB7 = cell(140,1);
phi = cell(140,1);
pd = cell(140,1);
phid = cell(140,1);
%Euler = zeros(3,140);

for i = t
    T07{i} = [RA{i},[AverageCurve(i,1);AverageCurve(i,2);AverageCurve(i,3)];0 0 0 1];
    TB7{i} = TB0*T07{i}; 
    phi{i} = EulerAngles(T07{i});
    pd{i} = [AverageCurve(i,1)'; AverageCurve(i,2)';AverageCurve(i,3)'];
    %Euler(:,i) = EulerRA{i};
end

disp('initial matrix');
disp(Forw_Kin(q(:,1)));
disp('goal matrix');
disp(T07{1});

%first point in 
pa = T07{1}(1:3,4); %wrt 0
xa = [pa;phi{1}]; 

disp('xa');
disp(xa);

K = .01;
xe_= zeros(length(xa), length(t));
e_ = zeros(length(xa), length(t));
ee = zeros(length(xa), length(t));
counter = 1;
jcounter = 0;

%xe(1:3) = position
%xe(4:6) = orientation

plim = 0.001;
olim = 0.0001;
tol = [plim; plim; plim; olim; olim; olim]; %error tolerance
for j = 1:3
    for i = 1:length(AverageCurve)
        
        xe = Forw_Kin(q(:,i));
        Ja = JacobianA(q(:,i));
        e = [pd{i};phi{i}] - xe;
      
        while (max(abs(e) > tol))
            Ja = JacobianA(qq);
            qdot = rem((pinv(Ja))*K*e,pi);
            qq = qq + qdot;
            xe = Forw_Kin(qq);

            e = [pd{i};phi{i}] - xe;
            e_(:,counter) = e;

            counter = counter + 1;
            disp(j);
            disp(i);
            disp(counter);
        end
        
        disp('step');
        disp(i);
        disp(j);
        jcounter = jcounter + 1;
        xe_(:,i) = xe;
        ee(:,i) = e;  
        qq = rem(qq,pi);
        q(:,i) = qq;
    end
end
disp('desired output of End Effector: ')
disp(xa)
disp('Final output of End Effector: ')
disp(xe_(:,length(xe_)))
disp('Computed Configuration: q')
disp(qq)
disp('Actual Transformation Matrix: (wrt 0)')
T_total = eye(4);

for i =1:7
    T = [cos(q(i,length(q))) -sin(q(i,length(q)))*cos(alpha(i)) sin(q(i,length(q)))*sin(alpha(i)) a(i)*cos(q(i,length(q)));...
         sin(q(i,length(q))) cos(q(i,length(q)))*cos(alpha(i)) -cos(q(i,length(q)))*sin(alpha(i)) a(i)*sin(q(i,length(q)));...
         0 sin(alpha(i)) cos(alpha(i)) d(i); 0 0 0 1];
     T_total = T_total*T;
end

disp(T_total)
dlmwrite('q.txt',q)


%plotting position 
figure
subplot(3,1,1);
plot(AverageCurve(:,1));
title("x: path");
subplot(3,1,2);
plot(AverageCurve(:,2));
title("y: path");
subplot(3,1,3);
plot(AverageCurve(:,3));
title("z: path");

figure
subplot(3,1,1);
plot((xe_(1,:)));
title("x: calculated");
subplot(3,1,2);
plot((xe_(2,:)));
title("y: calculated");
subplot(3,1,3);
plot((xe_(3,:)));
title("z: calculated");

%ORIENTATION 
figure
subplot(3,1,1);
plot((xe_(4,:)));
title("Phi(z)");
subplot(3,1,2);
plot((xe_(5,:)));
title("Theta(y)");
subplot(3,1,3);
plot((xe_(6,:)));
title("Psi(z)");

%plotting Joint Angles
figure
plot(q(1,:)*180/pi)
hold on
plot(q(2,:)*180/pi)
hold on
plot(q(3,:)*180/pi)
hold on
plot(q(4,:)*180/pi)
hold on
plot(q(5,:)*180/pi)
hold on
plot(q(6,:)*180/pi)
hold on
plot(q(7,:)*180/pi)

title("joint angles, deg")
legend('q1','q2','q3','q4','q5','q6','q7')

%plotting error at each frame
figure
plot(ee(1,:))
hold on
plot(ee(2,:))
hold on
plot(ee(3,:))
hold on
plot(ee(4,:))
hold on
plot(ee(5,:))
hold on
plot(ee(6,:))
hold on
plot(ones(size(ee(1,:)))*plim,'--k')
hold on
plot(ones(size(ee(1,:)))*plim*(-1),'--k')
hold on
plot(ones(size(ee(1,:)))*olim,'--k')
hold on
plot(ones(size(ee(1,:)))*olim*(-1),'--k')
hold on
title("xe error")
legend('e1','e2','e3','e4','e5','e6','plim','plim','olim','olim')

%plotting Angular Velocity of each joint
figure
plot(diff(q(1,:))*180/pi)
hold on
plot(diff(q(2,:))*180/pi)
hold on
plot(diff(q(3,:))*180/pi)
hold on
plot(diff(q(4,:))*180/pi)
hold on
plot(diff(q(5,:))*180/pi)
hold on
plot(diff(q(6,:))*180/pi)
hold on
plot(diff(q(7,:))*180/pi)
title("joint angles velocity, deg")
legend('q1','q2','q3','q4','q5','q6','q7')

