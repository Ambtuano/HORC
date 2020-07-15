function Euler = EulerAngles(T_total) %zyz
phi = atan2(-T_total(2,3),-T_total(1,3));
theta = atan2(-sqrt(T_total(1,3)^2 + T_total(2,3)^2),T_total(3,3));
psi = atan2(-T_total(3,2),T_total(3,1));

Euler = [phi; theta; psi];
 