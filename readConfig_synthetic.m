function config = readConfig_synthetic()
    config.algorithm = 'kdtree';
    config.trees = 8;
    config.checks = 64;
    
   
    config.datasetName = 'bunny10';
    config.flag = 3;       % 2或3，2表示利用距离，3表示利用面积
    config.maxIter = 100000;
    config.innerIter = 10;
    config.overlap=0.005;
%     config.overlap=0;
    config.k = 4;               % m in the paper
    config.pointPerSample = 8; % Nsample
    
%     config.epsilon = 0.15;
    config.epsilon = 0.005;
    config.nPairThresh = 2;   % Sample is accepted if the number of pairs that has distance > dDiffThreshold is greater than this
    config.nareaPairThresh =2;
    config.dDiffThresh = 0.5;
    config.dDiffThresh1 = 10; % Used for generating weight matrix
    config.pairDistThresh = 0.01;

         
    % For gsynthetic data   
    config.syntheticN = 5000;
    config.OutlierRate = 30;

    
    config.SamplesForCongruent = 500;
    config.extractPairEpsilon = 0.01;  % Grid size for extraing pairs
    
    config.angleDiffThreshold = 0.01; % For congruent testing
    
    %-------------DATASET-------------
    plyMaps = containers.Map(); plyMapsB = containers.Map();
    matMaps = containers.Map();
%     plyOutMaps = containers.Map(); 
    
    % Synthetic data with 10% outliers
%     plyMaps('bunny10') = './data/3DMatch/kitchen/cloud_bin_1.ply';
%     plyMapsB('bunny10') = './data/3DMatch/kitchen/cloud_bin_0.ply';
%     plyMaps('bunny10') = './data/ETH/gazebo_summer/Hokuyo_21.ply';
%     plyMapsB('bunny10') = './data/ETH/gazebo_summer/Hokuyo_0.ply';
%     plyMaps('bunny10') = './data/arch/s3.ply';
%     plyMapsB('bunny10') = './data/arch/s2.ply';
%     plyMaps('bunny10') = 'F:\新桌面\资料1\数字孪生\全新文献\点云\程序\数据\facade\s1.ply';
%     plyMapsB('bunny10') = 'F:\新桌面\资料1\数字孪生\全新文献\点云\程序\数据\facade\s2.ply';
%     plyMaps('bunny10') = 'F:\新桌面\资料1\数字孪生\全新文献\点云\程序\数据\bunny\data\bun045.ply';
%     plyMapsB('bunny10') = 'F:\新桌面\资料1\数字孪生\全新文献\点云\程序\数据\bunny\data\bun000.ply';
     plyMaps('bunny10') = 'F:\2021年\数字孪生\盘模型\pan3pcd.pcd';
    plyMapsB('bunny10') = 'F:\2021年\数字孪生\盘模型\real.ply';
%     plyMaps('bunny10') = 'C:\Users\23560\Desktop\data\Depth Long Throw\132554224282770898.ply';
%     plyMapsB('bunny10') = 'C:\Users\23560\Desktop\data\Depth Long Throw\132554224262770460.ply';
%     plyMaps('bunny10') = 'F:\新桌面\资料1\数字孪生\全新文献\点云\程序\数据\arch\s2.ply';
%     plyMapsB('bunny10') = 'F:\新桌面\资料1\数字孪生\全新文献\点云\程序\数据\arch\s4.ply';
    matMaps('bunny10') = './data/bunny10.mat';
   % plyOutMaps('bunny10') = '../Dataset/bunny/data/bunny10_synthetic';
    %%% More datasets can be addded in a similar manner
    
    config.plyPath = plyMaps(config.datasetName);
    config.plyPathB = plyMapsB(config.datasetName);
    config.matPath = matMaps(config.datasetName);
end