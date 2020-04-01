close all    % Close all open figures
clear        % Reset variables
clc          % Clear the command window

%LTEV2Vsim('help');

% Configuration file
configFile = 'BenchmarkPoisson.cfg';

% Simulation time (s)
T = 30;

% Beacon size (bytes)
B = 300;

%% LTE-V2V Controlled (BRAlgorithm 7)
% Network-controlled allocation Maximum Reuse Distance (MRD)
% Parameters:
% - Interval of scheduled reassignment: Treassign (s)

%LTEV2Vsim(configFile,'simulationTime',T,'BRAlgorithm',7,'Treassign',2,'Raw',150,...
%    'beaconSizeBytes',B);

%% LTE Autonomous (3GPP Mode 4)
% Autonomous allocation algorithm defined in 3GPP standard

LTEV2Vsim(configFile,'simulationTime',T,'BRAlgorithm',8,'Raw',150, ...
    'beaconSizeBytes',B, 'rho', 84, 'roadLength', 1800, 'maxPtx_dBm',23, 'minPtx_dBm',10, ...
    'MCS', 4);


