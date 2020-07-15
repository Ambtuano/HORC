%clear;
%close;
%Isolate the LTIB and LTIBA columns.
%Plot X, Y, Z over time and Y over Z for each one of the markers.


%importing marker position data from experiment
data = readtable('Subject_BH_2_27_2020_allysa01.csv');

%importing heel strike and coordinate data
frames = matfile('frames_of_heel_strikes_alyssa_exp_2_27_2020.mat');
heel_strikes = frames.frames_computed_heel_strikes;
lhee = matfile('lhee_xyz_coords_with_all_frames_alyssa_exp_2_27_2020.mat');
lhee_xyz = lhee.lhee_xyz_coordinates_with_all_frames;

%taking only LTIB and LTIBA data, 100hz 
t = str2double(data.Var1(3:end));

LTIBA_x = str2double(data.Subject_BH_2_27_2020_LTIBA(3:end));
LTIBA_y = str2double(data.Var25(3:end));
LTIBA_z = str2double(data.Var26(3:end));
LTIB_x = str2double(data.Subject_BH_2_27_2020_LTIB(3:end));
LTIB_y = str2double(data.Var28(3:end));
LTIB_z = str2double(data.Var29(3:end));
LHEE_x = lhee_xyz(:,2);
LHEE_y = lhee_xyz(:,3);
LHEE_z = lhee_xyz(:,4);
Position = [LTIBA_x LTIBA_y LTIBA_z LTIB_x LTIB_y LTIB_z LHEE_x LHEE_y LHEE_z];

%determine velocity by differentiating data
%multiply by 100 to get mm/s
LTIBA_dxdt = diff(Position(:,1),1,1)*100;
LTIBA_dydt = diff(Position(:,2),1,1)*100;
LTIBA_dzdt = diff(Position(:,3),1,1)*100;
LTIB_dxdt = diff(Position(:,4),1,1)*100;
LTIB_dydt = diff(Position(:,5),1,1)*100;
LTIB_dzdt = diff(Position(:,6),1,1)*100;
LHEE_dxdt = diff(Position(:,7),1,1)*100;
LHEE_dydt = diff(Position(:,8),1,1)*100;
LHEE_dzdt = diff(Position(:,9),1,1)*100;
Velocity = [LTIBA_dxdt LTIBA_dydt LTIBA_dzdt 
    LTIB_dxdt LTIB_dydt LTIB_dzdt 
    LHEE_dxdt LHEE_dydt LHEE_dzdt];

%I need to create an algorithm that decides whats the best heel strike
%frame to start at
%column 1 is position at each gait cycle
%column 2 is velocity at each gait cycle
gait_cycles = cell(134, 4);
for i = 1:length(heel_strikes(3:136))
    gait_cycles{i,1} = Position(heel_strikes(i+2):heel_strikes(i+3),:);
    gait_cycles{i,2} = Velocity(heel_strikes(i+2):heel_strikes(i+3),:);
end

%finding gait cycles with NaN values
%reminder 134 = length(vq) or length(gait_cycles)
rowNaN = zeros(length(gait_cycles),1);
colNaN = zeros(length(gait_cycles),1);

for i = 1:134
    row = find(isnan(gait_cycles{i,1}));
    rowNaN(i) = ~isempty(row);
end
%deleting cycles with NaN values
gait_cycles(logical(rowNaN),:) = [];


%vq is interpolating only position of gait_cycles to have the same length
%resample each cycle to 140 frames
%subtract the entire cycle of the first value of the cycle to prevent end
%effects in resampling
xq = 1:140;
vq = cell(length(gait_cycles),1);
for i = 1:1:size(Position,2)
    for j = 1:1:length(gait_cycles)
        y = resample(gait_cycles{j,1}(:,i)- gait_cycles{j,1}(1,i),140,length(gait_cycles{j,1}(:,i)));
        vq{j,i} = interp1(y,xq)' + gait_cycles{j,1}(1,i);
    end
end

%sum the position vectors to calculate average
% this will be used to hold the sum of all the points to calculate the average curve
vq2 = zeros(length(vq{1,1}), size(Position,2)); 
AverageCurve = zeros(length(vq{1,1}), size(Position,2)); 
for i = 1:1:size(Position,2)
    for j = 1:1:length(vq)
        vq2(:,i) = vq2(:,i) + vq{j,i}(:,1);
    end
    AverageCurve(:,i) = vq2(:,i)/length(vq{1});
end


%specific to this data 
%adjusts the coordinate system such that the robot will be on table and
%treadmill will be next to robot in the x direction (wrt global frame)
XCal = AverageCurve(1,1);
YCal = AverageCurve(1,2) + 500;
ZCal = AverageCurve(1,3);

for i = 1:3:size(AverageCurve,2)
    AverageCurve(:,i) = AverageCurve(:,i) - XCal; %x
    AverageCurve(:,i+1) = AverageCurve(:,i+1) - YCal; %y
    AverageCurve(:,i+2) = AverageCurve(:,i+2) -	ZCal; %z
end

% %Plots
%plot postion vs frame/over time for each individual step and then average]
Label = ["LTIBA_x (mm)", "LTIBA_y (mm)", "LTIBA_z (mm)", "LTIB_x (mm)", "LTIB_y (mm)", "LTIB_z (mm)", "LHEE_x (mm)", "LHEE_y (mm)", "LHEE_z (mm)"];
VarLabel = ["LTIBA", "LTIB", "LHEE"];
% 
% for x = 1:2
%     figure
%     for i = 1:size(Position,2)
%         for j = 1:1:length(gait_cycles)
%             subplot(size(Position,2)/2,2,i);
%             if x == 1
%                 plot(gait_cycles{j,1}(:,i));
%             end
%             if x == 2
%                 plot(AverageCurve(:,i)); 
%             end
%             hold on;
%             xlabel("# of frames");
%             ylabel(Label(i));
%             axis([0 130 0 1000]);
%             title('Position vs Time');
%         end
%     end
% end
% 
% % plot gait cycles of each market along yz axes
% figure
% subplot(3,1,1)
% for i = 1:size(Position,2)/3
% plot(AverageCurve(:,i*3-1),AverageCurve(:,i*3));
% hold on;
% xlabel("Y Position");
% ylabel("Z Position");
% title("Average Gait Cycle Curve: YZ Left view");
% axis([-750 0 -500 500]);
% end
% legend('LTIBA', 'LTIB', 'LHEE');
% 
% subplot(3,1,2)
% for i = 1:size(Position,2)/3
% plot(AverageCurve(:,i*3-2),AverageCurve(:,i*3-1));
% hold on;
% xlabel("X Position");
% ylabel("Y Position");
% title("Average Gait Cycle Curve: XY Top view");
% axis([-100 100 -1000 0]);
% end
% legend('LTIBA', 'LTIB', 'LHEE');
% 
% subplot(3,1,3)
% for i = 1:size(Position,2)/3
% plot(AverageCurve(:,i*3-2),AverageCurve(:,i*3));
% hold on;
% xlabel("X Position");
% ylabel("Z Position");
% title("Average Gait Cycle Curve: XZ Front view");
% axis([-100 100 -100 100]);
% end
% legend('LTIBA', 'LTIB', 'LHEE');

% figure
% plot3(AverageCurve(:,1),AverageCurve(:,2),AverageCurve(:,3),...
% AverageCurve(:,4),AverageCurve(:,5),AverageCurve(:,6),...
% AverageCurve(:,7),AverageCurve(:,8),AverageCurve(:,9));
% title("Average Gait Cycle Curve: XYZ");
% %axis([0 800 0 800 0 800]);
% legend('LTIBA', 'LTIB', 'LHEE');

%degrees (atan2(heel-tibia)) on YZ plane
thetaA_yz = (atan2(AverageCurve(:,8)-AverageCurve(:,2), AverageCurve(:,9)-AverageCurve(:,3)));%*180/pi);
thetaB_yz = (atan2(AverageCurve(:,8)-AverageCurve(:,5), AverageCurve(:,9)-AverageCurve(:,6)));%*180/pi);

%degrees (atan2(heel-tibia)) on XY 
thetaA_xy = (atan2(AverageCurve(:,7)-AverageCurve(:,1), AverageCurve(:,8)-AverageCurve(:,2)));%*180/pi);
thetaB_xy = (atan2(AverageCurve(:,7)-AverageCurve(:,4), AverageCurve(:,8)-AverageCurve(:,4)));%*180/pi);

%degrees (atan2(heel-tibia)) on XZ plane
thetaA_xz = (atan2(AverageCurve(:,7)-AverageCurve(:,1), AverageCurve(:,9)-AverageCurve(:,3)));%*180/pi);
thetaB_xz = (atan2(AverageCurve(:,7)-AverageCurve(:,4), AverageCurve(:,9)-AverageCurve(:,6)));%*180/pi);

%compensating for atan2 jump from -2pi to 0
for i = 1:length(AverageCurve)
    if thetaB_yz(i) < -pi/2
        thetaB_yz(i) = thetaB_yz(i) + 2*pi;
    end
end


%calculating Calibration Angle
thetaCalA_yz = mean(atan2(LHEE_y(1:700)-LTIBA_y(1:700),LHEE_z(1:700)-LTIBA_z(1:700)));%*180/pi;
thetaCalB_yz = mean(atan2(LHEE_y(1:700)-LTIB_y(1:700),LHEE_z(1:700)-LTIB_z(1:700)));%*180/pi;
thetaCalA_xy = mean(atan2(LHEE_x(1:700)-LTIBA_x(1:700),LHEE_y(1:700)-LTIBA_y(1:700)));%*180/pi;
thetaCalB_xy = mean(atan2(LHEE_x(1:700)-LTIB_x(1:700),LHEE_y(1:700)-LTIB_y(1:700)));%*180/pi;
thetaCalA_xz = mean(atan2(LHEE_x(1:700)-LTIBA_x(1:700),LHEE_z(1:700)-LTIBA_z(1:700)));%*180/pi;
thetaCalB_xz = mean(atan2(LHEE_x(1:700)-LTIB_x(1:700),LHEE_z(1:700)-LTIB_z(1:700)));%*180/pi;

phiX = thetaA_yz(:)-thetaCalA_yz;
thetaY = thetaA_xz(:)-thetaCalA_xz;
psiZ = thetaA_xy(:)-thetaCalA_xy;


figure

subplot(3,1,1);
plot((phiX)*180/pi); 
hold on;
plot((thetaB_yz-thetaCalB_yz)*180/pi);
xlabel("frame");
ylabel("angle (deg)");
title("Average Gait Cycle Curve: Calibrated YZ Plane, X axis" );
legend('LTIBA', 'LTIB');

subplot(3,1,2);
plot((thetaY)*180/pi); 
hold on;
plot((thetaB_xz-thetaCalB_xz)*180/pi);
xlabel("frame");
ylabel("angle (deg)");
title("Average Gait Cycle Curve: Calibrated XZ Plane, y axis" );
legend('LTIBA', 'LTIB');

subplot(3,1,3);
plot((psiZ)*180/pi); 
%hold on;
%plot((thetaB_xy-thetaCalB_xy)*180/pi);
xlabel("frame");
ylabel("angle (deg)");
title("Average Gait Cycle Curve: Calibrated XY Plane, z axis" );
legend('LTIBA');% 'LTIB');


%rotation matrix calculation, RPY(XYZ)
RA = cell(140,1);
EulerRA = cell(140,1);
for i = 1:length(RA)
    RA{i} = Rotmat(phiX(i),thetaY(i), psiZ(i));
    EulerRA{i} = EulerAngles(RA{i});
end
save('EulerRA.mat','EulerRA');
save('AverageCurve.mat','AverageCurve');
save('RA.mat','RA');
 