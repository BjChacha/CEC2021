%% 初始化
clear;
clc;

%% 定义参数
basePath = ".\IGDs";
baseAlgoIdx = 4;

%% 读取数据
srcList = dir(basePath);
algoNum = length(srcList);
i = 1;
maxTaskNum = 9;
while i <= algoNum
    if startsWith(srcList(i).name, ".") | startsWith(srcList(i).name, "_")
        srcList(i) = [];
        algoNum = algoNum - 1;
        continue
    end
    algoFolderList{i, 1} = srcList(i).folder + "\" + srcList(i).name;
    igdSrcList = dir(algoFolderList{i} + "\*.txt");
    igdSrcNum = length(igdSrcList);
    maxTaskNum = min(maxTaskNum, igdSrcNum);
    for j = 1:igdSrcNum
        x = getIGD(igdSrcList(j).folder + "\" + igdSrcList(j).name);
        igd{i,j} = x;
    end
    i = i + 1;
end

%% 统计
for i = 1:algoNum
    repeatedTimes = size(igd{i}, 2);
    for j = 1:maxTaskNum
        igd{i,j}(:,end+1) = min(igd{i,j}(:,1:repeatedTimes), [], 2);
        igd{i,j}(:,end+1) = max(igd{i,j}(:,1:repeatedTimes), [], 2);
        igd{i,j}(:,end+1) = median(igd{i,j}(:,1:repeatedTimes), 2);
        igd{i,j}(:,end+1) = std(igd{i,j}(:,1:repeatedTimes), 0, 2);
        igd{i,j}(:,end+1) = mean(igd{i,j}(:,1:repeatedTimes), 2);
    end
end

%% 质核检验
fprintf("基准算法： " + srcList(baseAlgoIdx).name + "\n");
fprintf("======== 质核检验 ========\n");
for i = 1:algoNum
    for j = 1:size(igd, 2)
        igd{i,j}(:, end+1) = 0;
    end
end
for i = 1:algoNum
    fprintf("---- " + srcList(i).name + " ----\n");
    for j = 1:size(igd, 2)
        if igd{i,j} == 0
            continue;
        end
        for k = 1:size(igd{i, j},1)
            [~, h] = ranksum(igd{baseAlgoIdx,j}(k,1:repeatedTimes), igd{i,j}(k,1:repeatedTimes));
            if h == 0
                igd{i,j}(k,end) = 0;
            elseif igd{baseAlgoIdx,j}(k,end-1) > igd{i,j}(k,end-1)
                igd{i,j}(k,end) = 1;
            else
                igd{i,j}(k,end) = -1;
            end
        end
        sim = sum(igd{i,j}(:,end) == 0);
        bet = sum(igd{i,j}(:,end) == 1);
        wor = sum(igd{i,j}(:,end) == -1);
        % fprintf('T%d ≈/-/+ : %d/%d/%d\n',j,sim,wor,bet); 
        fprintf('%d/%d/%d\n',sim,wor,bet);
    end
end

% 打印具体结果
fprintf("======== 具体结果 ========\n");
for i = 1:algoNum
    fprintf("---- " + srcList(i).name + " ----\n");
    for j = 1:size(igd, 2)
%         fprintf("T" + j + ": \n");
        for k = 1:size(igd{i, j},1)
            fprintf('%.2E',igd{i,j}(k,end-1));
            fprintf('±%.2E',igd{i,j}(k,end-2));
            switch igd{i,j}(k,end)
                case -1
                    fprintf('(-)');
                case 0
                    fprintf('(≈)');
                case 1
                    fprintf('(+)');
                otherwise
            end
            fprintf('\n');
        end 
%         fprintf('\n');
    end
end