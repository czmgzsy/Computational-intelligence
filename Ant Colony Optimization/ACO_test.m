%% 清空环境变量
clear all
clc

%% 导入数据
load citys_data.mat     % 从名为 'citys_data.mat' 的文件中加载数据

%31个城市的位置坐标(横坐标 纵坐标)保存在citys_data.mat文件中,变量citys为31行2列的数据
%第1列表示各个城市的横坐标,第2列表示各个城市的纵坐标

%% 计算城市间的相互距离
%城市间的距离矩阵是一个31阶方阵,由于对角线上的元素为0,而转移概率公式中的启发函数η是城市距离的倒数
%为了保证分母不为0,将对角线上的元素修正为一个非常小的正数(如10^-4或10^-5等)
n = size(citys,1);              % 获取城市数量
D = zeros(n,n);                 % 创建一个用于存储城市间距离的零矩阵
for i = 1:n
    for j = 1:n
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)) .^ 2));     %勾股定理
        else
            D(i,j) = 1 * 10^-4; % 对角线上的距离设置一个很小的值，避免出现零值
        end
    end
end

%% 初始化参数
m = 50;                         %蚂蚁数量
alpha = 1;                      %信息素重要程度因子
beta = 5;                       %启发函数重要程度因子
rho = 0.1;                      %信息素挥发因子
Q = 1;                          %蚁周模型中的常系数
Eta = 1 ./ D;                   %启发函数
Tau = ones(n,n);                %信息素矩阵,初始化为n行n列的全1方阵
Table = zeros(m,n);             %路径记录表,m只蚂蚁所以有m行记录,每只蚂蚁经过n个城市的次序所以是n列
iter = 1;                       %迭代次数初值
iter_max = 200;                 %最大迭代次数
Route_best = zeros(iter_max,n); %各代最佳路径
Length_best = zeros(iter_max,1);%各代最佳路径的长度
Length_ave = zeros(iter_max,1); %各代路径的平均长度


%% 迭代寻找最佳路径
while iter <= iter_max
    %随机产生各个蚂蚁的起点城市
    start = zeros(m,1);         % 创建一个初始城市编号的数组
    for i = 1:m
        temp = randperm(n);     %把1到n这些数随机打乱得到的一个数字序列
        start(i) = temp(1);     %每只蚂蚁的起始位置是随机的
    end
    Table(:,1) = start;         %路径记录表的第一列记为开始出发的城市
    
    %构建解空间
    citys_index = 1:n;          %生成了一个1至31的一维数组
    
    %逐个蚂蚁路径选择
    for i = 1 : m
        %逐个城市路径选择
        for j = 2 : n
            tabu = Table(i,1:(j - 1));                  %已访问的城市集合(禁忌表)
            allow_index = ~ismember(citys_index,tabu);  
            %找到尚未访问的城市,判断前一矩阵是否在后一矩阵中的逻辑值,并取反
            allow = citys_index(allow_index);           %待访问的城市集合
            P = allow;                                  % 创建概率数组，用于存储城市间转移概率
            %计算城市间转移概率
            for k = 1:length(allow)
                P(k) = Tau(tabu(end),allow(k)) ^ alpha * Eta(tabu(end),allow(k)) ^ beta;
            end
            P = P / sum(P);                             % 归一化概率
            %轮盘赌法选择下一个访问城市
            Pc = cumsum(P);
            target_index = find(Pc >= rand);            % 选择满足条件的城市编号
            target = allow(target_index(1));            % 获得下一个访问的城市编号
            Table(i,j) = target;                        % 将下一个城市编号放入路径记录表
        end
    end
    %计算各个蚂蚁的路径距离
    Length = zeros(m,1);                                % 创建数组，用于存储各个蚂蚁的路径距离
    for i = 1:m
        Route = Table(i,:);
        for j = 1:(n-1)
            Length(i) = Length(i) + D(Route(j),Route(j+1));% 计算路径长度
        end
        Length(i) = Length(i) + D(Route(n),Route(1));   % 考虑最后一个城市回到起点的距离
    end
    %计算最短路径距离及平均距离
    if iter == 1
        [min_Length,min_index] = min(Length);           % 获取最短路径长度及其索引(走出了最短路径的那只蚂蚁编号)
        Length_best(iter) = min_Length;                 % 存储本次迭代中的最短路径长度
        Length_ave(iter) = mean(Length);                % 计算本次迭代产生所有值的平均路径长度
        Route_best(iter,:) = Table(min_index);          % 存储最佳路径
    else
        [min_Length,min_index] = min(Length);           % 获取最短路径长度及其索引

        Length_best(iter) = min(Length_best(iter - 1),min_Length);
        % 更新最短路径长度,第一个参数是矩阵,第二个参数是标量,凡矩阵中的元素大于标量的,都用标量值替代

        Length_ave(iter) = mean(Length);                % 计算平均路径长度

        if Length_best(iter) == min_Length
            Route_best(iter,:) = Table(min_index,:);    % 存储最佳路径
        else
            Route_best(iter,:) = Route_best((iter-1),:);% 若未找到更优路径，保持上一代的最佳路径
        end
    end
    
    %更新信息素
    Delta_Tau = zeros(n,n); % 创建用于存储信息素增量的矩阵,n个城市所以是n阶
    
    %逐个蚂蚁计算
    for i = 1:m
        %逐个城市计算
        for j = 1:(n-1)
            Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q / Length(i);
            % 计算信息素增量
        end
        Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1)) + Q/Length(i);
        % 考虑回到起点的信息素增量
    end
    
    Tau = (1-rho) * Tau + Delta_Tau;    % 更新信息素
    
    %迭代次数+1,清空路径记录表
    iter = iter + 1;
    Table = zeros(m,n);
     
end

%% 结果显示
[Shortest_Length,index] = min(Length_best); % 获取最短路径长度及其索引
Shortest_Route = Route_best(index,:);         % 获取最短路径
disp(['最短距离:' num2str(Shortest_Length)]);% 显示最短路径长度
disp(['最短路径:' num2str([Shortest_Route Shortest_Route(1)])]);% 显示最短路径

%% 绘图
figure(1)
plot([citys(Shortest_Route,1);citys(Shortest_Route(1),1)],...
     [citys(Shortest_Route,2);citys(Shortest_Route(1),2)],'o-');
grid on
for i = 1:size(citys,1)
    text(citys(i,1),citys(i,2),['   ' num2str(i)]);
end
text(citys(Shortest_Route(1),1),citys(Shortest_Route(1),2),'       起点');
text(citys(Shortest_Route(end),1),citys(Shortest_Route(end),2),'       终点');
xlabel('城市位置横坐标')
ylabel('城市位置纵坐标')
title(['蚁群算法优化路径(最短距离:' num2str(Shortest_Length) ')'])
figure(2)
plot(1:iter_max,Length_best,'b',1:iter_max,Length_ave,'r:')
legend('最短距离','平均距离')
xlabel('迭代次数')
ylabel('距离')
title('各代最短距离与平均距离对比')




