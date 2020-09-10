%Calculates the end effector pose based on joint angles [7x1] using forward
%kinematics

function [ xe ] = Forw_Kin( q )
d = [.340 0 .400 0 .400 0 .126]';
alpha = [-pi/2 pi/2 pi/2 -pi/2 -pi/2 pi/2 0]';
a = zeros(7,1);

T_total = eye(4);
for i =1:7
    T = [cos(q(i)) -sin(q(i))*cos(alpha(i)) sin(q(i))*sin(alpha(i)) a(i)*cos(q(i));...
         sin(q(i)) cos(q(i))*cos(alpha(i)) -cos(q(i))*sin(alpha(i)) a(i)*sin(q(i));...
         0 sin(alpha(i)) cos(alpha(i)) d(i); 0 0 0 1];
     T_total = T_total*T;
end
xe = zeros(6,1);
xe(1:3) = T_total(1:3,4);
xe(4:6) = EulerAngles(T_total);
end

