% This script is written by Dengcheng Hao as a temoprary fix for the non-precise orbit 
% records that come with the data. It is written to delete the predicted orbit records
% that usually causes a abrrupt jump between themselves and the on-board collected records.
% This is not widely tested so use at your own risk!


% 此脚本解决LT-1A/B卫星原始xml文件中GPS数据丢失问题
% 需要输入参数包括filename_file_LED（使用GMTSAR提取的LED文件）
% use_number_point_to_hermit（对缺失数据插值时使用的数据量）
% outputFilename （hermit插值结果输出文件）

% 文件名
filename_file_LED = 'LT1.LED';
% 指定hermit插值所需数据量
use_number_point_to_hermit = 20;
% 指定输出文件名
outputFilename = 'output_LT1.txt';


%% 打开filename_file_LED文件
fileID = fopen(filename_file_LED, 'r');
if fileID == -1
    error('无法打开文件: %s', filename);
end
% 读取文件的前几行以获取元数据
headerLine = fgetl(fileID);
if isempty(headerLine)
    error('文件为空或格式不正确');
end
headerData = strsplit(headerLine);
numLines = str2double(headerData{1});
year = str2double(headerData{2});
day = str2double(headerData{3});
startSecond = str2double(headerData{4});
interval = str2double(headerData{5});
 
% 初始化存储数据的变量
data = [];
% 读取剩余的数据行
for i = 2:numLines+1 % +1 因为第一行是header
    line = fgetl(fileID);
    if isempty(line)
        error('文件中的数据行数少于预期');
    end
    lineData = strsplit(line);
    % 将字符串转换为数值
    yearData = str2double(lineData{1});
    dayData = str2double(lineData{2});
    secondData = str2double(lineData{3});
    x = str2double(lineData{4});
    y = str2double(lineData{5});
    z = str2double(lineData{6});
    vx = str2double(lineData{7});
    vy = str2double(lineData{8});
    vz = str2double(lineData{9});
    
    % 将数据添加到结构体数组中（或根据需要选择其他数据结构）
    newData = struct('Year', yearData, 'Day', dayData, 'Second', secondData, 'X', x, 'Y', y, 'Z', z, 'Vx', vx, 'Vy', vy, 'Vz', vz);
    data = [data; newData];
end
% 关闭文件
fclose(fileID);
% 将struct数据转为列表
data_x = transpose([data.X]);
data_vx = transpose([data.Vx]);
data_y = transpose([data.Y]);
data_vy = transpose([data.Vy]);
data_z = transpose([data.Z]);
data_vz = transpose([data.Vz]);

%% 判断Second 字段到一个数组中
seconds = [data.Second]; % 转换为列向量
sortedSeconds = sort(seconds);
 
diffs = diff(sortedSeconds);
 
% 定义一个阈值来检查连续性
tolerance = 1e-6; 
expectedInterval = 1; % 预期的时间间隔
 
% 检查差是否接近预期间隔，并记录不连续的索引（在排序后的数组中）
discontinuities = find(abs(diffs - expectedInterval) > tolerance) + 1;
 
% 输出结果
if isempty(discontinuities)
    disp('The Second values in data1 are continuous.');
else
    disp(['The Second values in data1 are not continuous. Discontinuities found at indices: ', num2str(discontinuities)]);
    
    
    % 如果需要，你可以输出不连续位置前后的 Second 值来帮助调试
    for idx = discontinuities
        %判断discontinuities是否位于头部，或者尾部
        i_small = idx - round((use_number_point_to_hermit + 1) / 2);
        if i_small <= 0
            i_small = 1;
        end
        i_max = idx + round((use_number_point_to_hermit + 1) / 2);
        if i_max > length(data)
            i_max = length(data);
        end
        % x插值
        hermit_second = transpose([data(i_small:i_max).Second]);
        hermit_x = transpose([data(i_small:i_max).X]);
        hermit_vx = transpose([data(i_small:i_max).Vx]);
        fprintf('Discontinuity at index %d: %f (prev) and %f (next)\n', idx, sortedSeconds(idx-1), sortedSeconds(idx));
        nodata = sortedSeconds(idx-1) + 1;
        f_x = Hermite(hermit_second,hermit_x,hermit_vx,nodata);
        fprintf('hermit x at index %d: %f \n', nodata, f_x);
        
        f_vx=interp1(hermit_second,hermit_vx,nodata,'spline');

        hermit_y = transpose([data(i_small:i_max).Y]);
        hermit_vy = transpose([data(i_small:i_max).Vy]);
        f_y = Hermite(hermit_second,hermit_y,hermit_vy,nodata);
        fprintf('hermit y at index %d: %f \n', nodata, f_y);

        f_vy=interp1(hermit_second,hermit_vy,nodata,'spline');

        hermit_z = transpose([data(i_small:i_max).Z]);
        hermit_vz = transpose([data(i_small:i_max).Vz]);
        f_z = Hermite(hermit_second,hermit_z,hermit_vz,nodata);
        fprintf('hermit z at index %d: %f \n', nodata, f_z);

        f_vz=interp1(hermit_second,hermit_vz,nodata,'spline');
        
        newStruct = struct('Year', yearData, 'Day', dayData, 'Second', nodata, 'X', double(f_x), 'Y', double(f_y), 'Z', double(f_z), 'Vx', f_vx, 'Vy', f_vy, 'Vz', f_vz);
        data1_cell = transpose(struct2cell(data));
        newEntry_cell = transpose(struct2cell(newStruct));
 

        data1_cell_expanded = [data1_cell(1:idx-1, :); newEntry_cell; data1_cell(idx:length(data), :)];
 
        data_fields = fieldnames(data(1)); 
        data1_expanded = cell2struct(data1_cell_expanded, data_fields, 2);
        data = data1_expanded;
    end

    % 结果输出
    % 打开文件以写入
    fileID = fopen(outputFilename, 'w');
    % 检查文件是否成功打开
    if fileID == -1
        error('无法打开文件: %s', outputFilename);
    end
    numRows = height(data1_cell_expanded);
    numCols = width(data1_cell_expanded);
    
    rowStr_01 = '';
    rowStr_01 = [rowStr_01,sprintf('%.0f', numRows),' ',sprintf('%.0f', data1_cell_expanded{1,1}),' ',sprintf('%.0f', data1_cell_expanded{1,2}), ' ',sprintf('%.3f', data1_cell_expanded{1,3}), ' ',sprintf('%.6f', 1),' '];
    fprintf(fileID,'%s\n', rowStr_01);
    for i = 1:numRows
        rowStr = '';
        rowStr = [rowStr, sprintf('%.0f', data1_cell_expanded{i,1}),' ',sprintf('%.0f', data1_cell_expanded{i,2}), ' ',sprintf('%.3f', data1_cell_expanded{i,3}),' ',sprintf('%.6f', data1_cell_expanded{i,4}),' ',sprintf('%.6f', data1_cell_expanded{i,5}),' ',sprintf('%.6f', data1_cell_expanded{i,6}),' ',sprintf('%.6f', data1_cell_expanded{i,7}),' ',sprintf('%.6f', data1_cell_expanded{i,8}),' ',sprintf('%.6f', data1_cell_expanded{i,9}),' '];
    
        fprintf(fileID,'%s\n', rowStr);
    end
     
    % 关闭文件（如果使用 fopen 打开了文件）
    fclose(fileID);
end


function f = Hermite( x,y,y_1,x0 )
%Hermite插值函数
%   x为已知数据点的x坐标
%   y为已知数据点的y坐标
%   y_1为数据点y值导数
%   x0为插值点的x坐标
syms t;
f = 0.0;
if(length(x) == length(y))
    if(length(x) == length(y_1))
        n = length(x);
    else
        disp('y和y的导数维数不相等');
        renturn;
    end
else
    disp('x和y的维数不相等');
    return;
end
%以上为输入判断和确定“n”的值
for i = 1:n
    h =  1.0;
    a =  0.0;
    for j = 1:n
        if(j ~= i)
            h = h*(t-x(j))^2/((x(i)-x(j))^2);%求得值为(li(x))^2
            a = a + 1/(x(i)-x(j));   %求得ai（x）表达式之中的累加部分
        end
    end
    
    f = f + h*((x(i)-t)*(2*a*y(i)-y_1(i))+y(i));
    
    if(i == n)
        if(nargin == 4)
            f = subs(f, 't' , x0);  %输出结果
        else
            f = vpa(f,6);  %输出精度为有效数字为6位的函数表达式
        end
    end
end
end

