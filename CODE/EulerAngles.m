%Calculates Euler ZYZ angles based on the given transformation matrix [4x4]

function Euler = EulerAngles(T_total) %zyz
phi = atan2(-T_total(2,3),-T_total(1,3));
theta = atan2(-sqrt(T_total(1,3)^2 + T_total(2,3)^2),T_total(3,3));
psi = atan2(-T_total(3,2),T_total(3,1));

if psi > 0
    psi = psi - 2*pi;
end
Euler = [phi; theta; psi];
 