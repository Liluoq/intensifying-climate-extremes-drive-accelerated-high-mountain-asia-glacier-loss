%% CryoSat-2数据处理详细说明文档
%
% 参考：revised_ICESat2_code/Processing_Documentation.m
% 本脚本提供CryoSat-2数据处理的技术细节说明

%% 1. 核心技术差异：11×11邻域平均

fprintf('=== CryoSat-2关键技术特征 ===\n\n');

% 邻域尺寸对比
neighbor_comparison = [
    "数据源", "足迹尺寸", "邻域大小", "像元数", "代码";
    "ICESat-2", "~17m", "3×3", "9", "for j=-1:1; for k=-1:1";
    "CryoSat-2", "~300m", "11×11", "121", "for j=-5:5; for k=-5:5"
];

disp('邻域尺寸对比:');
disp(neighbor_comparison);

fprintf('\n计算过程:\n');
fprintf('CryoSat-2足迹: ~300m\n');
fprintf('NASADEM分辨率: 30m\n');
fprintf('所需像元数: 300m/30m ≈ 10\n');
fprintf('邻域尺寸: 11×11（确保覆盖）\n');
fprintf('总像元数: 11×11 = 121\n\n');

%% 2. 代码实现示例

fprintf('=== 11×11邻域平均代码实现 ===\n\n');

% 模拟示例
fprintf('示例代码:\n');
fprintf('-----------------------------------------------------\n');
fprintf('neighbor_size = 11;\n');
fprintf('half_size = 5;\n');
fprintf('\n');
fprintf('for i = 1:size(data, 1)\n');
fprintf('    [row, col] = latlon2pix(ref_srtm, data(i,2), data(i,3));\n');
fprintf('    row = ceil(row);\n');
fprintf('    col = ceil(col);\n');
fprintf('    \n');
fprintf('    ele = 0.0;\n');
fprintf('    for j = -5:5  %% 11个像元\n');
fprintf('        for k = -5:5  %% 11个像元\n');
fprintf('            ele = ele + double(srtm(row+j, col+k));\n');
fprintf('        end\n');
fprintf('    end\n');
fprintf('    \n');
fprintf('    data(i, 23) = ele / 121;  %% 平均\n');
fprintf('end\n');
fprintf('-----------------------------------------------------\n\n');

%% 3. 数据格式说明

fprintf('=== 数据格式说明 ===\n\n');

% 第一阶段输出格式
stage1_format = [
    "列号", "字段", "说明";
    "1", "年份", "从文件名提取";
    "2", "纬度", "WGS84";
    "3", "经度", "WGS84";
    "4", "观测高程", "CryoSat-2测高值";
    "5", "时间", "相对2000-01-01的天数";
    "6", "大地水准面", "Geoid高度";
    "7", "地表高程", "elevation - geoid";
    "8", "HMA22区域", "1-22";
    "9", "GTN区域", "RGI区域";
    "10", "0.5°网格", "网格编号";
    "11", "1°网格", "网格编号";
    "12", "TP标识", "青藏高原";
    "13", "冰川标识", "0=冰川，1=非冰川"
];

disp('第一阶段输出格式（HMA_CryoSat2_glacier.txt）:');
disp(stage1_format);

% 第二阶段输出格式
stage2_format = [
    "列号", "字段", "说明";
    "23", "NASADEM高程", "11×11邻域平均（关键！）";
    "24", "高程变化", "地表高程 - NASADEM(11×11)"
];

disp('第二阶段新增列（HMA_CryoSat2_update.txt）:');
disp(stage2_format);

%% 4. 质量控制策略

fprintf('\n=== 质量控制策略 ===\n\n');

quality_control_steps = [
    "步骤", "过滤条件", "阈值";
    "1. 高程范围", "观测高程", "<8900m";
    "2. 边界过滤", "HMA边界", "inpolygon";
    "3. 冰川掩膜", "RGI掩膜", "glacier_flag=0";
    "4. 高程变化", "abs(dh)", "<400m";
    "5. 百分位数", "dh范围", "5%-85%";
    "6. 中位数过滤", "abs(dh-median)", "<75m"
];

disp('质量控制步骤:');
disp(quality_control_steps);

%% 5. 时间处理

fprintf('\n=== 时间处理说明 ===\n\n');

fprintf('时间参考基准: 2000-01-01\n');
fprintf('数据时段: 2010-2021\n');
fprintf('时间单位: 天\n\n');

fprintf('时间转换公式:\n');
fprintf('  NetCDF时间戳（秒） → 天数:\n');
fprintf('  days = floor(time_20_ku / 86400)\n\n');
fprintf('  天数 → 日期:\n');
fprintf('  date = datetime(2000,1,1) + days(days_value)\n\n');

fprintf('年度时间窗口:\n');
for year = 2010:2:2020
    year_start = datenum(year, 1, 1) - datenum(2000, 1, 1);
    year_end = datenum(year+1, 1, 1) - datenum(2000, 1, 1);
    fprintf('  %d年: 天数 %d - %d\n', year, year_start, year_end);
end

%% 6. 统计方法

fprintf('\n=== 鲁棒统计方法 ===\n\n');

fprintf('中位数过滤法:\n');
fprintf('  med = median(dh);\n');
fprintf('  valid = abs(dh - med) < 75;\n');
fprintf('  mean_robust = mean(dh(valid));\n\n');

fprintf('百分位数法:\n');
fprintf('  lower = prctile(dh, 5);\n');
fprintf('  upper = prctile(dh, 85);\n');
fprintf('  valid = (dh >= lower) & (dh <= upper);\n');
fprintf('  保留约88%%的数据\n\n');

%% 7. 与ICESat-2的完整对比

fprintf('=== CryoSat-2 vs ICESat-2 完整对比 ===\n\n');

full_comparison = [
    "特征", "CryoSat-2", "ICESat-2";
    "测高原理", "雷达高度计", "激光高度计";
    "足迹尺寸", "~300m", "~17m";
    "NASADEM邻域", "11×11 (121)", "3×3 (9)";
    "数据时段", "2010-2021", "2018-2024+";
    "时间参考", "2000-01-01", "2018-01-01";
    "轨道重复", "369天", "91天";
    "穿透深度校正", "不需要", "需要";
    "数据密度", "较低", "较高";
    "质量控制", "相对宽松", "相对严格"
];

disp(full_comparison);

%% 8. 处理流程图

fprintf('\n=== 数据处理流程 ===\n\n');

fprintf('L2I/TEMPO NetCDF 文件\n');
fprintf('        ↓\n');
fprintf('[extract_cryosat2_tracks] - 提取RGI冰川区域\n');
fprintf('        ↓\n');
fprintf('HMA_CryoSat2_glacier.txt\n');
fprintf('        ↓\n');
fprintf('[enhance_cryosat2_data] - 11×11邻域NASADEM\n');
fprintf('        ↓\n');
fprintf('HMA_CryoSat2_update.txt\n');
fprintf('        ↓\n');
fprintf('[compute_cryosat2_statistics]\n');
fprintf('        ↓\n');
fprintf('分析结果（MAT）\n');
fprintf('        ↓\n');
fprintf('[generate_cryosat2_plots]\n');
fprintf('        ↓\n');
fprintf('图表（PNG/FIG）\n\n');

%% 9. 常见问题

fprintf('=== 常见问题解答 ===\n\n');

fprintf('Q1: 为什么CryoSat-2用11×11邻域？\n');
fprintf('A1: 因为CryoSat-2足迹~300m，远大于ICESat-2的~17m，\n');
fprintf('    需要更大邻域来匹配其空间分辨率。\n\n');

fprintf('Q2: 能否用3×3邻域处理CryoSat-2？\n');
fprintf('A2: 不能！会导致空间不匹配，结果不可靠。\n\n');

fprintf('Q3: NetCDF读取失败怎么办？\n');
fprintf('A3: 检查变量名是否正确（L2I用_poca_20_ku后缀）。\n\n');

fprintf('Q4: 如何验证邻域平均是否正确？\n');
fprintf('A4: 检查列23的NASADEM值应该是平滑的，\n');
fprintf('    且与单点NASADEM值接近但更平滑。\n\n');

%% 10. 完成提示

fprintf('=== 文档说明完成 ===\n');
fprintf('请运行主脚本: CryoSat2_Main_Processing_Pipeline.m\n');
fprintf('关键记住: 11×11邻域，不是3×3！\n');
