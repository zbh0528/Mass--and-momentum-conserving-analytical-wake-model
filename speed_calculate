function Actual_wind_speed = speed_calculate(trans_xy_position, wf, turbine, ind_v)
% 计算风电场中每台风机的实际风速，考虑尾流效应
% 输入参数:
%   trans_xy_position - 风机的坐标位置矩阵 (2 x turbine_num)
%   wf                - 风场参数结构体，包含风机和环境参数
%   turbine           - 涡轮机参数结构体，包含风机数量、涡轮半径、湍流强度等
%   ind_v             - 当前风速索引
%
% 输出参数：
%   Actual_wind_speed - 各台风机的实际风速向量
[~, sorted_index] = sort(-trans_xy_position(2,:));  % 按 y 坐标逆序排列
T = turbine.turbine_num;  % 获取风机的数量
% 初始化每台风机的湍流强度和实际风速
I = ones(T, 1) * turbine.turbulence_intensity;  % 初始湍流强度
Actual_wind_speed = ones(T, 1) * wf.velocity(ind_v);  % 初始实际风速

% 遍历每台风机，计算其有效风速，考虑尾流效应
for j = 2:T
    idxJ = sorted_index(j);  % 第 j 台风机的排序索引
    delta_I = 0;  % 初始化湍流系数增长量
    TjSpeed = Actual_wind_speed(idxJ);  % 当前风机的初始风速
    
    % 考虑之前的风机对当前风机的尾流影响
    for i = 1:j-1
        idxI = sorted_index(i);  % 前一台风机的排序索引
        x_dist = trans_xy_position(1, idxJ) - trans_xy_position(1, idxI);  % 风机间 x 方向距离
        y_dist = trans_xy_position(2, idxJ) - trans_xy_position(2, idxI);  % 风机间 y 方向距离
        
        % 判断风机 i 是否位于风机 j 之前，并且风速处于涡轮机工作的范围内
        if y_dist < 0  && TjSpeed >= turbine.CutIn && TjSpeed <= turbine.CutOut
            % 通过样条插值计算风机 j 的风力系数 Ct
            TiSpeed = Actual_wind_speed(idxI);  % 当前风机的初始风速
            Ct = spline(turbine.data.ws, turbine.data.Ct, TiSpeed);
            % Ct = 0.7;
            % 计算尾流潜在核心长度 x0
            x0 = (1 + sqrt(1 - Ct)) / (sqrt(2) * (4 * 0.58 * I(idxI) + 2 * 0.077 * sqrt(1 - Ct))); 
            
            % 计算尾流增长率 k
            k = 0.003678 + 0.3837 * I(idxI);
            
            % 计算尾流宽度 sigma
            sigma = k * (abs(y_dist) - x0) + turbine.rotor_diameter / sqrt(8);
            
            % 计算尾流导致的风速损失 delta_U
            delta_U = TjSpeed * (1 - sqrt(1 - Ct / (8 * (sigma / turbine.rotor_diameter)^2))) ...
                * exp(-((x_dist^2 + (-turbine.hub_height)^2) / (2 * sigma^2)));
            
            % 更新风机 j 的风速
            TjSpeed = TjSpeed - delta_U;
            
            % 计算湍流强度增长量 delta_I
            R = turbine.rotor_radius + k * abs(y_dist);  % 尾流半径
            Aw = cal_interaction_area(abs(x_dist), abs(y_dist), turbine.rotor_radius, R);  % 交互面积
            delta_ij = (0.73 * ((1 - sqrt(1 - Ct)) / 2)^0.8325 * I(idxI)^-0.0325 * ...
                (abs(y_dist) / turbine.rotor_diameter)^-0.32) * (Aw * 4 / (pi * turbine.rotor_diameter^2));
            
            % 选择最大的 delta_I 作为湍流强度变化量
            if delta_I < delta_ij
                delta_I = delta_ij;
            end
        end
    end
    
    % 更新当前风机的实际风速
    Actual_wind_speed(idxJ) = TjSpeed;
    
    % 更新湍流强度
    I(idxJ) = sqrt(delta_I^2 + turbine.turbulence_intensity^2);
end
end
