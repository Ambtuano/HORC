
AverageCurvemat = load('AverageCurve.mat');
AverageCurve = AverageCurvemat.AverageCurve; %mm

% figure 
% TibA = animatedline; 
% TibB = animatedline;
% Lhee = animatedline;
% axis([0 800 0 800]); ]
% for i = 1:length(AverageCurve)
%     addpoints(TibA,AverageCurve(i,2),AverageCurve(i,3)) 
%     addpoints(TibB,AverageCurve(i,5),AverageCurve(i,6)) 
%     addpoints(Lhee,AverageCurve(i,8),AverageCurve(i,9))
%     xlabel("y position (mm)");
%     ylabel("z position (mm)");
%     title("Average Gait Cycle Curve: Left side perspective" );
%     drawnow
% end
% hold on

% plot gait cycles of each marketr along each axes
% figure
% subplot(3,1,1)
% for i = 1:size(Position,2)/3
% plot(AverageCurve(:,i*3-2));
% hold on;
% title("Average Gait Cycle Curve: x view");
% axis([0 140 -100 100]);
% end
% legend('LTIBA', 'LTIB', 'LHEE');
% 
% subplot(3,1,2)
% for i = 1:size(Position,2)/3
% plot(AverageCurve(:,i*3-1));
% hold on;
% title("Average Gait Cycle Curve: Y view");
% axis([0 140 -1000 100]);
% end
% legend('LTIBA', 'LTIB', 'LHEE');
% 
% subplot(3,1,3)
% for i = 1:size(Position,2)/3
% plot(AverageCurve(:,i*3));
% hold on;
% title("Average Gait Cycle Curve: Z Front view");
% axis([0 140 -500 100]);
% end
% legend('LTIBA', 'LTIB', 'LHEE');




% %2d plot xy
% for i = 1:length(AverageCurve)
%     plot(AverageCurve(:,1),AverageCurve(:,2),...
%     AverageCurve(:,4),AverageCurve(:,5),...
%     AverageCurve(:,7),AverageCurve(:,8));
%     hold on
%     plot(AverageCurve(i,1),AverageCurve(i,2),'*',...
%     AverageCurve(i,4),AverageCurve(i,5),'*',...
%     AverageCurve(i,7),AverageCurve(i,8),'*');
% 
%     title("Average Gait Cycle Curve: XY Top");
%     axis([0 800 0 800 0 800]);
%     legend('LTIBA', 'LTIB', 'LHEE');
%     pause(.001);
%     hold off;
% end


%2d plot yz
for j = 1:4
    for i = 1:length(AverageCurve)
        plot(AverageCurve(:,2),AverageCurve(:,3),...
        AverageCurve(:,5),AverageCurve(:,6),...
        AverageCurve(:,8),AverageCurve(:,9));
        hold on
        plot(AverageCurve(i,2),AverageCurve(i,3),'*',...
        AverageCurve(i,5),AverageCurve(i,6),'*',...
        AverageCurve(i,8),AverageCurve(i,9),'*');

        title("Average Gait Cycle Curve: YZ left view");
        axis([-600 100 -600 100]);
        legend('LTIBA', 'LTIB', 'LHEE');
        disp(i);
        pause(.001);
        hold off;
    end
end

% %2d plot xz
% for i = 1:length(AverageCurve)
%     plot(AverageCurve(:,1),AverageCurve(:,3),...
%     AverageCurve(:,4),AverageCurve(:,6),...
%     AverageCurve(:,7),AverageCurve(:,9));
%     hold on
%     plot(AverageCurve(i,1),AverageCurve(i,3),'*',...
%     AverageCurve(i,4),AverageCurve(i,6),'*',...
%     AverageCurve(i,7),AverageCurve(i,9),'*');
% 
%     title("Average Gait Cycle Curve: XZ front");
%     axis([0 800 0 800 0 800]);
%     legend('LTIBA', 'LTIB', 'LHEE');
%     pause(.001);
%     hold off;
% end

% %3d plot
% figure
%     plot3(AverageCurve(:,1),AverageCurve(:,2),AverageCurve(:,3),'*',...
%     AverageCurve(:,4),AverageCurve(:,5),AverageCurve(:,6),'*',...
%     AverageCurve(:,7),AverageCurve(:,8),AverageCurve(:,9),'*',...
%     0,0,0,'o');
%     title("Average Gait Cycle Curve: XYZ");
%     legend('LTIBA', 'LTIB', 'LHEE');
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     