%% 《计算智能》——张军 清华大学出版社  教材中第六章对PSO算法思想来源的介绍节选
% Wilson 在20世纪70年代就指出:"至少在理论上,在群体觅食的过程中,群体中的每一个个体
% 都会受益于所有成员在这个过程中所发现和积累的经验." 因此PSO直接采用了这一思想
% PSO算法在决策过程中使用两类重要的信息:一是自身的经验,二是其他人的经验
%% 清空命令行窗口和工作区
clc
clear

%% 目标函数 求解此函数在[0,50]区间上的最小值
f = @(x) x.*sin(x).*cos(2*x) - 2*x.*sin(3*x)+3*x.*sin(4*x);
% @是函数句柄的符号,这句话的意思是说 f 是关于 x 的函数,后边是函数的表达式
% help function_handle

%% 初始化算法参数
population = 50;            % 粒子的数量
dimension = 1;              % 解的维度
iter_max = 100;             % 迭代终止条件
location_limit = [0,50];    % 粒子运动区间限制,本例中即为自变量区间限制
velocity_limit = [-10,10];  % 速度限制
omega = 0.7;                % 粒子惯性权重
c1 = 0.5;                   % 自身经验系数
c2 = 0.5;                   % 群体经验系数

%% 绘制目标函数图像
figure(1)
ezplot(f,[location_limit(1),location_limit(2)]);
% 绘制目标函数的图像,区间由location_limit(1)(2)控制

%% 初始化粒子群
pop = location_limit(1)+(location_limit(2)-location_limit(1))*rand(population,dimension);
% 粒子群的初始位置是随机生成的,分布在定义域上
% pop是一个population行,dimension列的矩阵

v = rand(population,dimension);                 % 随机生成粒子群的初始速度
optimal_position1 = pop;                        % 粒子个体的历史最佳位置
optimal_position2 = zeros(1,dimension);         % 粒子群体的历史最佳位置
optimal_value1 = ones(population,dimension)*inf;% 粒子个体的历史最佳函数值
optimal_value2 = ones(1,dimension)*inf;         % 粒子群体的历史最佳函数值

optimal_record = ones(iter_max,1) * inf;        % 记录每轮迭代中群体最佳函数值

hold on
plot(pop,f(pop),'ro');
title('初始状态图');
hold off

pause(0.5);

figure(2);

%% 迭代循环
for iter = 1:iter_max
    % 1.函数值计算及最优记录更新
    f_value = f(pop);                       % 粒子当前函数值
    for i = 1 : population
        if f_value(i) < optimal_value1(i)   % 粒子当前函数值小于粒子历史最佳函数值
            optimal_value1(i) = f_value(i); % 更新粒子历史最佳函数值
            optimal_position1(i,:)=pop(i,:);% 更新粒子历史最佳位置
        end  
    end
    if min(optimal_value1) < optimal_value2
        [optimal_value2,index]= min(optimal_value1);% 得到粒子群最优函数值及其索引
        optimal_position2 = optimal_position1(index,:);
        % 粒子的历史最佳位置优于群体的历史最佳位置时,将群体的历史最佳位置更新
    end
    
    % 2.更新粒子的速度和位置
    v = v*omega ...                                         % 初速度*惯性权重
    + c1*rand*(optimal_position1-pop) ...                   % 自身经验
    + c2*rand*(repmat(optimal_position2,population,1)-pop); % 群体经验
    % 公式比较长,用三个点,表示续行换行命令
    % 语法有B = repmat(A,m,n)，即把 A 作为 B 的元素，B 由 m×n 个 A 平铺而成
    % optimal_position2记录的是群体最优值,其中仅有一条记录
    % 要对整个群体做矩阵运算,需要先将optimal_position2扩充为与整个群体相同的维度
    
    % 3.速度越界处理
    v(v<velocity_limit(1)) = velocity_limit(1);
    v(v>velocity_limit(2)) = velocity_limit(2);

    % 4.更新粒子位置
    pop = pop + v;

    % 5.位置越界处理
    pop(pop<location_limit(1)) = location_limit(1);
    pop(pop>location_limit(2)) = location_limit(2);
    
    % 6.迭代过程中最优值记录
    optimal_record(iter) = optimal_value2;      % 记录本轮迭代中的群体最优值
    
    % 7.绘图
    x = location_limit(1):0.01:location_limit(2);

    subplot(1,2,1)  % 将对话窗分割为1行2列,并对画幅的第1列进行操作
    plot(x,f(x),'b-',pop,f(pop),'ro');

    subplot(1,2,2)
    plot(optimal_record);
    title('最优适应度进化过程');
    pause(0.01);
    iter = iter + 1;
    
end


x = location_limit(1):0.01:location_limit(2);
figure(3);
plot(x,f(x),'b-',pop,f(pop),'ro');
title('最终位置状态');
disp('最优函数值',num2str(optimal_value2));
disp('\n最优解',num2str(optimal_position2));



