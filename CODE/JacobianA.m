%Calculates analytical jacobian [6x6] through given joint angles (q)

function [ Ja ] = JacobianA( q )
d = [.340 0 .400 0 .400 0 .126]';
alpha = [-pi/2 pi/2 pi/2 -pi/2 -pi/2 pi/2 0]';
a = zeros(7,1);

T = cell(7,1);
for i = 1:7
    T{i} = [cos(q(i)) -sin(q(i))*cos(alpha(i)) sin(q(i))*sin(alpha(i)) a(i)*cos(q(i));...
            sin(q(i)) cos(q(i))*cos(alpha(i)) -cos(q(i))*sin(alpha(i)) a(i)*sin(q(i));...
            0 sin(alpha(i)) cos(alpha(i)) d(i); 0 0 0 1];
end

J = zeros(6,7);
T01 = T{1};
T02 = T01*T{2};
T03 = T02*T{3};
T04 = T03*T{4};
T05 = T04*T{5};
T06 = T05*T{6};
T07 = T06*T{7};

J(1:3,1) = cross([0;0;1], (T07(1:3,4) - [0;0;0]));
J(1:3,2) = cross(T01(1:3,3), (T07(1:3,4)-T01(1:3,4)));
J(1:3,3) = cross(T02(1:3,3), (T07(1:3,4)-T02(1:3,4)));
J(1:3,4) = cross(T03(1:3,3), (T07(1:3,4)-T03(1:3,4)));
J(1:3,5) = cross(T04(1:3,3), (T07(1:3,4)-T04(1:3,4)));
J(1:3,6) = cross(T05(1:3,3), (T07(1:3,4)-T05(1:3,4)));
J(1:3,7) = cross(T06(1:3,3), (T07(1:3,4)-T06(1:3,4)));

J(4:6,1) = [0;0;1];
J(4:6,2) = T01(1:3,3);
J(4:6,3) = T02(1:3,3);
J(4:6,4) = T03(1:3,3);
J(4:6,5) = T04(1:3,3);
J(4:6,6) = T05(1:3,3);
J(4:6,7) = T06(1:3,3);

xe = Forw_Kin(q);
phi = xe(4);
theta = xe(5);
Ta = [1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 0 -sin(phi) cos(phi)*sin(theta);
    0 0 0 0 cos(phi) sin(phi)*sin(theta);
    0 0 0 1 0 cos(theta)];
Ja = pinv(Ta)*J; %analytical Jacobian


end