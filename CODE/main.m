clear all;
close all;

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


%initial angular configuration
q0 = zeros(7,1);
qq = q0;
q = zeros(length(q0), length(t));
q(:,1) = q0;

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

% d = [.340 0 .400 0 .400 0 .126]'; %m
% alpha = [-pi/2 pi/2 pi/2 -pi/2 -pi/2 pi/2 0]';
% a = zeros(7,1);
% 
% for i =1:7
%     T{i} = [cos(qq(i)) -sin(qq(i))*cos(alpha(i)) sin(qq(i))*sin(alpha(i)) a(i)*cos(qq(i));
%          sin(qq(i)) cos(qq(i))*cos(alpha(i)) -cos(qq(i))*sin(alpha(i)) a(i)*sin(qq(i));
%          0 sin(alpha(i)) cos(alpha(i)) d(i); 0 0 0 1];
% end
% 
TB0 = [0 0 1 0; 1 0 0 0; 0 1 0 0; 0 0 0 1];
% J = zeros(6,7);
% T01 = T{1}; 
% T02 = T01*T{2}; 
% T03 = T02*T{3};
% T04 = T03*T{4};
% T05 = T04*T{5};
% T06 = T05*T{6};
% T07 = T06*T{7};
T07 = [RA{1},[AverageCurve(1,1);AverageCurve(1,2);AverageCurve(1,3)];0 0 0 1];
TB7 = TB0*T07; 

disp('initial matrix');
disp(Forw_Kin(q(:,1)));
disp('goal matrix');
disp(T07);

phi = EulerAngles(T07);

%first point in 
pa = T07(1:3,4); %wrt 0
xa = [pa;phi]; 

disp('xa');
disp(xa);

pd = [AverageCurve(t,1)'; AverageCurve(t,2)'];
phid = AverageCurve(t,3)';
K = .1;
xe_= zeros(length(xa), length(t));
e_ = zeros(length(xa), length(t));
counter = 1;
for i = 1:length(t)
    xe = Forw_Kin(q(:,i));
    Ja = JacobianA(q(:,i));
    e = [pd(:,i); phid(i);phi] - xe;
    disp(i);
    disp(e);
    while (max(abs(e)) > 0.00001)
        Ja = JacobianA(qq);
        qdot = (pinv(Ja))*K*e;
        qq = qq + qdot;
        
        xe = Forw_Kin(qq);
        e = [pd(:,i); phid(i);phi] - xe;
        e_(:,counter) = e;
        
        counter = counter + 1;
        disp(counter);
        disp(i);
    end
    xe_(:,i) = xe;
    q(:,i+1) = qq;
end

disp('desired output of End Effector: ')
disp(xa)
disp('Final output of End Effector: ')
disp(xe)
disp('Computed Configuration: ')
disp(qq)
disp('Actual Transformation Matrix: ')
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

disp(T_total)
dlmwrite('q.txt',q)


