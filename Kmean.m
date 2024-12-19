% K-Means 聚类算法动画演示
% 作者：ChatGPT
% 日期：2024-04-27

clear; clc; close all;

%% 参数设置
numPoints = 300;    % 数据点数量
K = 3;               % 簇的数量
maxIterations = 10;  % 最大迭代次数
rng(42);             % 设置随机种子以便结果可重复

%% 生成随机数据
% 生成 K 个不同高斯分布的数据簇
data = [];
centers = [5, 5; 15, 15; 25, 5]; % 三个簇的中心
for i = 1:K
    data = [data; mvnrnd(centers(i,:), [3, 3], numPoints/K)];
end

% 绘制初始数据
figure('Color', 'w');
scatter(data(:,1), data(:,2), 30, 'k', 'filled');
title('K-均值聚类动画演示');
xlabel('X 轴');
ylabel('Y 轴');
hold on;

%% 初始化质心
% 随机选择 K 个数据点作为初始质心
randomIndices = randperm(size(data,1), K);
centroids = data(randomIndices, :);
% 绘制初始质心
hCentroids = scatter(centroids(:,1), centroids(:,2), 100, 'r', 'x', 'LineWidth', 2);
legend('数据点', '初始质心');

%% K-Means 迭代
for iter = 1:maxIterations
    % 分配步骤：将每个点分配到最近的质心
    distances = pdist2(data, centroids, 'euclidean');
    [~, labels] = min(distances, [], 2);
    
    % 更新步骤：重新计算质心
    newCentroids = zeros(K, size(data,2));
    for j = 1:K
        if sum(labels == j) == 0
            % 如果某个簇没有分配到任何点，随机重新初始化
            newCentroids(j, :) = data(randi(size(data,1)), :);
        else
            newCentroids(j, :) = mean(data(labels == j, :), 1);
        end
    end
    
    % 绘制当前聚类结果
    clf; % 清空当前图形
    hold on;
    
    colors = lines(K);
    for j = 1:K
        clusterData = data(labels == j, :);
        scatter(clusterData(:,1), clusterData(:,2), 30, colors(j,:), 'filled');
    end
    % 绘制质心
    hCentroids = scatter(newCentroids(:,1), newCentroids(:,2), 100, 'r', 'x', 'LineWidth', 2);
    
    % 图形设置
    title(sprintf('K-均值聚类动画演示 - 迭代 %d', iter));
    xlabel('X 轴');
    ylabel('Y 轴');
    legendEntries = arrayfun(@(x) sprintf('簇 %d', x), 1:K, 'UniformOutput', false);
    legend([hCentroids; scatter(nan, nan)], ['质心', legendEntries], 'Location', 'best');
    grid on;
    axis([0 30 0 20]);
    drawnow;
    
    % 检查收敛性
    if isequal(centroids, newCentroids)
        disp(['算法在迭代 ' num2str(iter) ' 时收敛。']);
        break;
    end
    
    centroids = newCentroids;
    pause(1); % 暂停一秒以观察动画
end

%% 最终结果
% 最终聚类结果已经在最后一次迭代中显示
disp('K-均值聚类完成。');

%% 进一步示例：确定最佳K值并使用K-Means++初始化

% 以下代码展示了如何使用肘部法和轮廓系数来确定最佳K值，并使用K-Means++初始化方法

% MATLAB 示例：确定最佳K值并使用K-Means++初始化

clear; clc; close all;

%% 生成随机数据
rng(42); % 设置随机种子
numPoints = 300;
K_true = 4; % 实际簇数（用于生成数据）
data = [];
centers_true = [5, 5; 15, 15; 25, 5; 35, 15];
for i = 1:K_true
    data = [data; mvnrnd(centers_true(i, :), [3, 3], numPoints/K_true)];
end

%% 肘部法确定K值
wcss = [];
K_range = 1:10;
for K = K_range
    [~, ~, sumd] = kmeans(data, K, 'Replicates', 10, 'Display', 'off');
    wcss(K) = sum(sumd);
end

% 绘制肘部法图
figure('Color', 'w');
plot(K_range, wcss, '-o', 'LineWidth', 2);
title('肘部法确定最佳K值');
xlabel('簇的数量 K');
ylabel('簇内平方和 (WCSS)');
grid on;

%% 计算轮廓系数
silhouette_avg = zeros(1, length(K_range));
for K = K_range
    if K ==1
        silhouette_avg(K) = 0;
    else
        labels = kmeans(data, K, 'Replicates', 10, 'Display', 'off');
        silhouette_avg(K) = mean(silhouette(data, labels));
    end
end

% 绘制轮廓系数图
figure('Color', 'w');
plot(K_range, silhouette_avg, '-s', 'LineWidth', 2);
title('轮廓系数确定最佳K值');
xlabel('簇的数量 K');
ylabel('平均轮廓系数');
grid on;

%% 选择K值（例如，根据肘部法和轮廓系数，选择K=4）
K_optimal = 4;

%% 使用K-Means++初始化
centroids_init = kmeans_pp_init(data, K_optimal);

%% 运行K-均值算法
opts = statset('Display','final','MaxIter',300);
[idx, centroids, sumd, D] = kmeans(data, K_optimal, 'Start', centroids_init, ...
    'Replicates', 1, 'Options', opts); % 将 'Replicates' 设置为 1

%% 可视化聚类结果
figure('Color', 'w');
colors = lines(K_optimal);
hold on;
for k = 1:K_optimal
    cluster_data = data(idx == k, :);
    scatter(cluster_data(:,1), cluster_data(:,2), 30, colors(k,:), 'filled');
end
scatter(centroids(:,1), centroids(:,2), 100, 'kx', 'LineWidth', 2);
title(sprintf('K-均值聚类结果 (K=%d)', K_optimal));
xlabel('X 轴');
ylabel('Y 轴');
legend_entries = arrayfun(@(x) sprintf('簇 %d', x), 1:K_optimal, 'UniformOutput', false);
legend([legend_entries, {'质心'}], 'Location', 'best');
grid on;
hold off;

%% K-Means++ 初始化函数
function centroids = kmeans_pp_init(data, K)
    N = size(data, 1);
    centroids = zeros(K, size(data, 2));
    % 随机选择第一个质心
    idx = randi(N);
    centroids(1, :) = data(idx, :);
    
    for k = 2:K
        distances = min(pdist2(data, centroids(1:k-1, :), 'euclidean'), [], 2);
        probabilities = distances.^2 / sum(distances.^2);
        cumulative_prob = cumsum(probabilities);
        r = rand();
        idx = find(cumulative_prob >= r, 1, 'first');
        centroids(k, :) = data(idx, :);
    end
end
