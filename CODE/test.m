clear all
close all
clc

vrep = remApi('remoteApi');
vrep.simxFinish(-1);
clientID = vrep.simxStart('127.0.0.1', 19999, true, true, 5000, 5);

%insert q
q = load('q.txt','q');
if (clientID>-1)
    disp('Connected to remote API server')
    
    h = [0,0,0,0,0,0,0];
    
    [r,h(1)] = vrep.simxGetObjectHandle(clientID, 'LBR_iiwa_7_R800_joint1', vrep.simx_opmode_blocking);
    [r,h(2)] = vrep.simxGetObjectHandle(clientID, 'LBR_iiwa_7_R800_joint2', vrep.simx_opmode_blocking);    
    [r,h(3)] = vrep.simxGetObjectHandle(clientID, 'LBR_iiwa_7_R800_joint3', vrep.simx_opmode_blocking);
    [r,h(4)] = vrep.simxGetObjectHandle(clientID, 'LBR_iiwa_7_R800_joint4', vrep.simx_opmode_blocking);    
    [r,h(5)] = vrep.simxGetObjectHandle(clientID, 'LBR_iiwa_7_R800_joint5', vrep.simx_opmode_blocking);
    [r,h(6)] = vrep.simxGetObjectHandle(clientID, 'LBR_iiwa_7_R800_joint6', vrep.simx_opmode_blocking);    
    [r,h(7)] = vrep.simxGetObjectHandle(clientID, 'LBR_iiwa_7_R800_joint7', vrep.simx_opmode_blocking);
    
    %while true
    for i = 1:7
        vrep.simxSetJointTargetPosition(clientID, h(i), q(i,1), vrep.simx_opmode_streaming);
    end
    pause(10);
        for j = 1:length(q)
            for i = 1:7
                vrep.simxSetJointTargetPosition(clientID, h(i), q(i,j), vrep.simx_opmode_streaming);
            end
            pause(0.001);
        end
    for i = 1:7
        vrep.simxSetJointTargetPosition(clientID, h(i), q(i,length(q)), vrep.simx_opmode_streaming);
    end    
    %end
    
else
       disp('Failed connection to remote API server');
end
vrep.delete();
disp('program ended');