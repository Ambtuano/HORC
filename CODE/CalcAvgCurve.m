%Attempt at creating generalized code of CalcAverageCurve.m just using
%input data


function [RA, AverageCurve] = CalcAvgCurve( LTIBA, LHEE, heel_strikes ) % Assumes coordinate and frame information of markers

% Information Expected 
%{
Format: LHEE and LTIBA in similar formats, .mat files with columns to
be [frames, x coordinate, y coordinate, z coordinate] respectively 

lhee = matfile('lhee_xyz_coords.mat'); 
ltiba = matfile('ltiba_xyz_coords.mat'); 
frames = matfile('frames_of_heel_strikes.mat'); 

LTIBA = ltiba.lhee_xyz_coordinates_with_all_frames; 
LHEE = lhee.lhee_xyz_coordinates_with_all_frames; 
heel_strikes = frames.frames_computed_heel_strikes;
%}

%{

lines for testing:

lhee = matfile('lhee_xyz_coords_with_all_frames_alyssa_exp_2_27_2020.mat');
LHEE = lhee.lhee_xyz_coordinates_with_all_frames;
LTIBA_x = str2double(data.Subject_BH_2_27_2020_LTIBA(3:end));
LTIBA_y = str2double(data.Var25(3:end));
LTIBA_z = str2double(data.Var26(3:end));
LTIBA = [zeros(length(LTIBA_x),1), LTIBA_x, LTIBA_y,LTIBA_z]
frames = matfile('frames_of_heel_strikes_alyssa_exp_2_27_2020.mat');
heel_strikes = frames.frames_computed_heel_strikes;
[RA, AverageCurve] = CalcAvgCurve( LTIBA, LHEE, heel_strikes );

%}
%taking only LTIB and LTIBA data, 100hz 
AvgHS = round(mean(diff(heel_strikes)));

LTIBA_x = LTIBA(:,2);
LTIBA_y = LTIBA(:,3);
LTIBA_z = LTIBA(:,4);
LHEE_x = LHEE(:,2);
LHEE_y = LHEE(:,3);
LHEE_z = LHEE(:,4);
Position = [LTIBA_x LTIBA_y LTIBA_z LHEE_x LHEE_y LHEE_z];

%determine velocity by differentiating data
%multiply by 100 to get mm/s
LTIBA_dxdt = diff(Position(:,1),1,1)*100;
LTIBA_dydt = diff(Position(:,2),1,1)*100;
LTIBA_dzdt = diff(Position(:,3),1,1)*100;
LHEE_dxdt = diff(Position(:,4),1,1)*100;
LHEE_dydt = diff(Position(:,5),1,1)*100;
LHEE_dzdt = diff(Position(:,6),1,1)*100;
Velocity = [LTIBA_dxdt LTIBA_dydt LTIBA_dzdt LHEE_dxdt LHEE_dydt LHEE_dzdt];

%I need to create an algorithm that decides whats the best heel strike
%frame to start at
%column 1 is position at each gait cycle
%column 2 is velocity at each gait cycle
gait_cycles = cell(134, 1);
for i = 1:length(heel_strikes(3:136))
    gait_cycles{i,1} = Position(heel_strikes(i+2):heel_strikes(i+3),:);
    gait_cycles{i,2} = Velocity(heel_strikes(i+2):heel_strikes(i+3),:);
end

%finding gait cycles with NaN values
%reminder 134 = length(vq) or length(gait_cycles)
rowNaN = zeros(length(gait_cycles),1);

for i = 1:134
    row = find(isnan(gait_cycles{i,1}));
    rowNaN(i) = ~isempty(row);
end
%deleting cycles with NaN values
gait_cycles(logical(rowNaN),:) = [];


%vq is interpolating only position of gait_cycles to have the same length
%resample each cycle to Average Gait Cycle number of frames
%subtract the entire cycle of the first value of the cycle to prevent end
%effects in resampling
xq = 1:AvgHS;
vq = cell(length(gait_cycles),1);
for i = 1:1:size(Position,2)
    for j = 1:1:length(gait_cycles)
        y = resample(gait_cycles{j,1}(:,i)- gait_cycles{j,1}(1,i),AvgHS,length(gait_cycles{j,1}(:,i)));
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
%Label = ["LTIBA_x (mm)", "LTIBA_y (mm)", "LTIBA_z (mm)", "LTIB_x (mm)", "LTIB_y (mm)", "LTIB_z (mm)", "LHEE_x (mm)", "LHEE_y (mm)", "LHEE_z (mm)"];
%VarLabel = ["LTIBA", "LTIB", "LHEE"];
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
% 3D Plot of average curves
% figure
% plot3(AverageCurve(:,1),AverageCurve(:,2),AverageCurve(:,3),...
% AverageCurve(:,4),AverageCurve(:,5),AverageCurve(:,6));
% title("Average Gait Cycle Curve: XYZ");
% %axis([0 800 0 800 0 800]);
% legend('LTIBA', 'LTIB', 'LHEE');

%degrees (atan2(heel-tibia (z,y))) on YZ plane
theta_yz = (atan2( AverageCurve(:,6)-AverageCurve(:,3),AverageCurve(:,5)-AverageCurve(:,2)));%*180/pi);

%degrees (atan2(heel-tibia (x,y))) on XY 
theta_xy = mean(atan2(AverageCurve(:,4)-AverageCurve(:,1), AverageCurve(:,5)-AverageCurve(:,2)));%*180/pi);

%degrees (atan2(heel-tibia (x,z))) on XZ plane
theta_xz = mean(atan2(AverageCurve(:,4)-AverageCurve(:,1), AverageCurve(:,6)-AverageCurve(:,3)));%*180/pi);


%compensating for atan2 jump from -2pi to 0
for i = 1:length(AverageCurve)
    if theta_yz(i) < 0
        theta_yz(i) = theta_yz(i) + pi;
    end
end

% calculating Calibration Angle based on standing posture
thetaCalA_yz = mean(atan2(LHEE_z(1:700)-LTIBA_z(1:700),LHEE_y(1:700)-LTIBA_y(1:700)));%*180/pi;
thetaCalA_xy = mean(atan2(LHEE_x(1:700)-LTIBA_x(1:700),LHEE_y(1:700)-LTIBA_y(1:700)));%*180/pi;
thetaCalA_xz = mean(atan2(LHEE_x(1:700)-LTIBA_x(1:700),LHEE_z(1:700)-LTIBA_z(1:700)));%*180/pi;


phiZ = theta_xy(:)-thetaCalA_xy;
thetaY = theta_xz(:)-thetaCalA_xz;
psiX = theta_yz(:)-thetaCalA_yz;

figure
plot((psiX)*180/pi); 
xlabel("frame");
ylabel("angle (deg)");
title("Average Gait Cycle Curve: Calibrated YZ Plane, X axis" );
legend('LTIBA', 'LTIB');

%rotation matrix calculation, RPY(XYZ)
RA = cell(AvgHS,1);
EulerRA = cell(AvgHS,1);
for i = 1:length(RA)
    RA{i} = Rotmat(phiZ, thetaY, psiX(i));
    EulerRA{i} = EulerAngles(RA{i});
end

save('AverageCurve.mat','AverageCurve'); %Position
save('RA.mat','RA'); %Rotation Matrix
