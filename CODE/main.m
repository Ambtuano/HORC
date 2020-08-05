clear all;
close all;
%last updated: 8.4.2020

%given information
%frequency = 100hz
%Gait Cycle = 140 frames = 1.4s
t = 1:1:140;

%approximate link diameter:
width = .068*2; %m

AverageCurvemat = load('AverageCurve.mat');
AverageCurve = AverageCurvemat.AverageCurve/1000; %m
RAmat = load('RA.mat');
RA = RAmat.RA;
EulerRAmat = load('EulerRA.mat');
EulerRA = EulerRAmat.EulerRA;
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
Euler = zeros(3,140);

for i = t
    T07{i} = [RA{i},[AverageCurve(i,1);AverageCurve(i,2);AverageCurve(i,3)];0 0 0 1];
    TB7{i} = TB0*T07{i}; 
    phi{i} = EulerAngles(T07{i});
    pd{i} = [AverageCurve(i,1)'; AverageCurve(i,2)'];
    phid{i} = AverageCurve(i,3)';
    Euler(:,i) = EulerRA{i};
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

for j = 1:2
    for i = 1:length(AverageCurve)
        xe = Forw_Kin(q(:,i));
        Ja = JacobianA(q(:,i));
        e = [pd{i}; phid{i};phi{i}] - xe;
        while (max(abs(e(1:3))) > 0.0001 && max(abs(e(4:6))) > 0.001)
            Ja = JacobianA(qq);
            qdot = rem((pinv(Ja))*K*e,pi);
            qq = qq + qdot;
            %check qq within -pi,pi -> normalize
            xe = Forw_Kin(qq);

            e = [pd{i}; phid{i};phi{i}] - xe;
            e_(:,counter) = e;

            counter = counter + 1;
            disp(counter);
            disp(i);
            disp(j);
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
qq = rem(qq,pi);
disp('desired output of End Effector: ')
disp(xa)
disp('Final output of End Effector: ')
disp(xe_(:,length(xe_)))
disp('Computed Configuration: q')
disp(qq)
disp('Actual Transformation Matrix: (wrt 0)')
d = [.3105 0 .400 0 .390 0 .083]';
alpha = [pi/2 -pi/2 -pi/2 pi/2 pi/2 -pi/2 0]';
a = zeros(7,1);
T_total = eye(4);

for i =1:7
    T = [cos(qq(i)) -sin(qq(i))*cos(alpha(i)) sin(qq(i))*sin(alpha(i)) a(i)*cos(qq(i));
         sin(qq(i)) cos(qq(i))*cos(alpha(i)) -cos(qq(i))*sin(alpha(i)) a(i)*sin(qq(i));
         0 sin(alpha(i)) cos(alpha(i)) d(i); 0 0 0 1];
     T_total = T_total*T;
end

disp(TB0*T_total)
dlmwrite('q.txt',q)


%position 

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
plot(q(1,:))
hold on
plot(q(2,:))
hold on
plot(q(3,:))
hold on
plot(q(4,:))
hold on
plot(q(5,:))
hold on
plot(q(6,:))
hold on
plot(q(7,:))

title("joint angles, rad")
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

title("xe error")
legend('e1','e2','e3','e4','e5','e6')



