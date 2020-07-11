clear all;
close all;

%Isolate the LTIB and LTIBA columns.
%Plot X, Y, Z over time and Y over Z for each one of the markers.

%importing marker position data from experiment
data = readtable('Subject_BH_2_27_2020_allysa01.csv');

%importing heel strike and coordinate data
frames = matfile('frames_of_heel_strikes_alyssa_exp_2_27_2020.mat');
heel_strikes = frames.frames_computed_heel_strikes;
lhee = matfile('lhee_xyz_coords_with_all_frames_alyssa_exp_2_27_2020.mat');
lhee_xyz = lhee.lhee_xyz_coordinates_with_all_frames;

%taking only LTIB and LTIBA data
t = str2double(data.Var1(3:end));
LTIBA_y = str2double(data.Var25(3:end));
LTIB_y = str2double(data.Var28(3:end));
Position = [LTIBA_y LTIB_y];

%determine velocity by differentiating data
%multiply by 100 to get mm/s
LTIBA_dydt = diff(Position(:,1),1,1)*100;
LTIB_dydt = diff(Position(:,2),1,1)*100;

%removing outliers
A = isoutlier(LTIB_dydt);
LTIBA_dydt(A)=[];
LTIBA_y(A)=[];
B = isoutlier(LTIB_dydt);
LTIB_dydt(B)=[];
LTIB_y(B)=[];
Velocity = [ LTIBA_dydt LTIB_dydt];
Position = [LTIBA_y LTIB_y];

% %Truncating data pt1
% %determining data from treadmill velocity < 50 m/s (v20)
% %using heel strike data points to get entire gait cycles
% %v20 ends at heel_strikes(3) as after the third heel strike we can see the
% %gait cycle curve become larger indicating an increase in the subjects
% %walking speed and therefore an increase in treadmill speed
% v20 = zeros(heel_strikes(3),2);
% for i = 1:heel_strikes(3)
%     v20(i,2) = Velocity(i,2);
%     v20(i,1) = Velocity(i,1);
% end

% %Truncating data pt2
% %determining data from treadmill velocity = 50 m/s 
% %by removing v20 datapoints
% v50 = zeros(length(Velocity(heel_strikes(3):heel_strikes(137),1)),2);
% for i = 1:1:length(v50)
%     v50(i,1) = Velocity(heel_strikes(3)+i-1,1);
%     v50(i,2) = Velocity(heel_strikes(3)+i-1,2);
% end

%splitting data
%split the velocity data such that we have a single curve per gait cycle
%each row is an individual gait cycle.
%each column represents position(1) or velocity(2)
%within each cell, there is nx2 double, where each column is LTIBA and LTIB
%respectively and each row is a singular frame from that gait cycle
%the first 2 gait cycles are excluded as this was the point in time where
%the treadmill velocity was below 50m/s


gait_cycles = cell(134, 2);

for i = 1:length(heel_strikes(3:136))
    gait_cycles{i,1} = Position(heel_strikes(i+2):heel_strikes(i+3),:);
    gait_cycles{i,2} = Velocity(heel_strikes(i+2):heel_strikes(i+3),:);
end

Label_x2 = ["LTIBA_y (mm)", "LTIB_y (mm)"];
Label_y2 = ["LTIBA_d_y_d_t (mm/s)", "LTIB_d_y_d_t (mm/s)"];

% %plot phase plots
% figure 
% for i = 1:1:length(Label_x2)
%     for j = 1:1:length(heel_strikes(3:136))
%         subplot(1,2,i);
%         plot(gait_cycles{j,1}(:,i),gait_cycles{j,2}(:,i));
%         hold on;
%         ylabel(Label_y2(i));
%         xlabel(Label_x2(i));
%         
%         title('Phase plots');
%     end
% end

xq = 1:130;
vq = cell(134,1);
for i = 1:1:length(Label_x2)
    for j = 1:1:length(heel_strikes(3:136))
        vq{j,i} = interp1(gait_cycles{j,1}(:,i),xq)';
    end
end

%plot postion vs frame
figure 
for i = 1:2
    for j = 1:1:length(heel_strikes(3:136))
        subplot(1,2,i);
        plot(gait_cycles{j,1}(:,i));
        hold on;
        xlabel("# of frames");
        ylabel(Label_x2(i));
        axis([0 130 200 700]);
        title('position over time curve');
    end
end



%sum the position vectors to calculate average
% this will be used to hold the sum of all the points to calculate the average curve
vq2 = zeros(130,2); 

for i = 1:1:length(Label_x2)
    for j = 1:1:length(vq{1,1})
        vq2(:,i) = vq2(:,i) + vq{j,i}(:,1);
    end
end

AverageCurve(:,1) = vq2(:,1)/length(vq{1,1});
AverageCurve(:,2) = vq2(:,2)/length(vq{1,1});

figure 
for i = 1:1:length(Label_x2)
    for j = 1:1:length(heel_strikes(3:136))
        subplot(1,2,i);
        plot(xq,AverageCurve(:,i));
        xlabel("# of frames");
        ylabel(Label_x2(i));
        axis([0 130 200 700]);
        title('position over time average curve');
    end
end



