%% 清空工作区
clc
clear all;

%% 绘制目标函数
% y = 10 * sin(5 * x)+7 * abs(x - 5)+10, x=[0,10]
figure(1);
ezplot('10*sin(5*x)+7*abs(x-5)+10',[0,10]);
hold on
xlabel('自变量/x');
ylabel('应变量/y');

%% GA算法参数定义
popsize = 20;               % 种群大小
chromlength = 20;           % 二进制编码长度
p_cross = 0.7;              % 交叉概率
p_mutation = 0.1;          % 变异概率
iter_max = 500;             % 最大迭代次数
bestfit = -inf;             % 初始化适应度值为-∞
bestindividual = 0;         % 初始化最优解
result = zeros(2,iter_max); % 对存储每轮迭代最优解的result矩阵预分配内存,提高运行速度
%% 初始化种群
pop = round(rand(popsize,chromlength));
% rand函数用0~1之间的随机数生成了一个popsize行,chromlength列的矩阵
% round函数将此矩阵四舍五入为0或1的整数

%% 循环迭代
for iter = 1:iter_max
    % 1.首先将二进制转为十进制
    [px,py] = size(pop);                    % 种群的行列数
    for j = 1:py
        pop1(:,j) = pop(:,j) * 2^(py-j);    % 按位权值展开
    end
    temp = sum(pop1,2);                     % 行求和sum(*,2)
    pop2 = 0+temp*((10-0)/2^chromlength);

    % 2.计算函数目标值
    fitvalue = 10*sin(5*pop2)+7*abs(pop2-5)+10;

    % 3.选择操作(轮盘赌法)
    p_fitvalue = fitvalue / sum(fitvalue);  % 单个适应度与所有适应度之和的比
    p_fitvalue = cumsum(p_fitvalue);        % 离散的概率值转为0~1区间,累计和
    
    p_rand = sort(rand(px,1));              % px行,1列个由小到大排序的随机数
    % 用p_rand去与做完了累计和的p_fitvalue做比较
    % 适应度越高就有越大的几率在p_fitvalue的累计和中占据更大的区间
    % 将会有较多的p_rand值落在其间
    % 用这种方法保留下了适应度值较大的基因
    i = 1;      % 指向随机数的指针
    j = 1;      % 指向累积和区间的指针
    while i <= px
        if p_rand(i) < p_fitvalue(j)
            newpop(i,:) = pop(j,:);    %适应度高的基因以较大概率保留到了新种群中
            i = i + 1;
        else
            j = j + 1;  %小于当前选择概率区间的随机数"用完了",才让概率区间后移一位
        end
    end
    
    % 4.交叉操作(单点交叉)
    % newpop种群已经经过了选择操作,现对其进行交叉操作
    % 并将交叉后的种群保存进newpop2矩阵中
    for i = 1:2:px-1 % 让i与i+1进行交叉,所以i遍历染色体中编号为奇数的就可以了
        if rand < p_cross
            point = round(rand * py);   % 在py位二进制数中随机选一个交叉点位
            newpop2(i,:) = [newpop(i,1:point),newpop(i+1,point+1:end)];
            newpop2(i+1,:) = [newpop(i+1,1:point),newpop(i,point+1:end)];
        else
            newpop2(i,:) = newpop(i,:);
            newpop2(i+1,:) = newpop(i+1,:);
            % 一定几率不发生交叉
        end
    end
    
    % 5.变异操作
    for i = 1:px
        if rand < p_mutation                    % 参数定义时的变异概率
            mutation_point = round(rand*py);   % 发生变异时的点位
            if mutation_point <= 0
                mutation_point = 1;             % 检测越界
            end
            newpop3(i,:) = newpop2(i,:);        % 变异后的种群存储在newpop3矩阵中
            if newpop3(i,mutation_point) == 0   % 将变异点位上的二进制数取反
                newpop3(i,mutation_point) = 1;
            else
                newpop3(i,mutation_point) = 0;
            end
        else                                    % p_mutation设置为一个很小的数
            newpop3(i,:) = newpop2(i,:);        % 大部分newpop2都没有发生变异
        end
    end

    % 6.更新种群,重新计算种群的适应度,更新Best记录
    pop = newpop3;
    % 和第一步时一样,先将二进制转换为十进制,计算出小分划区间中的实数值
    [px,py] = size(pop);
    for j = 1:py
        pop1(:,j) = 2^(py-j) * pop(:,j); % 二进制数按位权值展开
    end  % chromlength位二进制数能表示的十进制数范围为[0,2^chromlength-1]
    temp = sum(pop1,2);     % 按位权展开的十进制数行求和
    pop2 = 0 + ((10 - 0)/2^chromlength) * temp; % 将小分划的编号映射为定义域里的实数值

    fitvalue = 10*sin(5*pop2)+7*abs(pop2-5)+10; 
    % 将pop2矩阵代入目标函数计算适应度

    for i = 1:px
        if bestfit < fitvalue(i)
            bestfit = fitvalue(i);      % 找到本轮迭代的最优f(x)
            bestindividual = pop2(i);   % 找到本轮迭代的最优x
        end
    end

    result(1,iter) = bestindividual;    % result矩阵第1行存储每轮迭代的最优x
    result(2,iter) = bestfit;           % result矩阵第2行存储每轮迭代的最优f(x)
    plot(result(1,iter),result(2,iter),'b*','MarkerSize',5)
    title(strcat('当前迭代代数:',num2str(iter)));
    pause(0.03);
end

plot(result(1,end),result(2,end),'r*','MarkerSize',10);
grid on
hold off
fprintf(['The best x is ---> \t',num2str(bestindividual), ...
    '\nThe fitness value is --->',num2str(bestfit),'\n']);





