% 原有代码保持不变...
%% ---------------- Clear All --------------------------------------------
clear all
close all
clc
clear

%% 鍑芥暟璇存槑
% TD_window 涓烘埅宄板苟鎷熷悎閫艰繎鐢紱
% time_pad 涓烘硶澶ц幓鑷达細閫夋嫨

%% 材料和数据参数
c0 = 0.000299792458 * 10^6;     % [um*THz] 光速
Z0 = 376.73031;               % [V/A] 自由空间的阻抗
eps0 = 1 / (Z0 * c0);         % [A/(V*um*THz)] 自由空间的介电常数
smp_thickness = 1 * 1e3;    % [um] 样品厚度
n1 = 1;
nn = 1;                    % 基板的折射率

windowing_width = 5;
type = 'gauss';  % 时间窗口类型
Zero_padding = 2^13; % 零填充长度
Num_eff_0 = 10;
Num_eff_1 = 1230; 
%% 读取Excel文件
excelPath = 'D:\张雨林\导电力学协同优化\原始数据\太赫兹吸收测试\修正后TFSI1mm202508281118.xls';
output_folder = 'D:\张雨林\导电力学协同优化\原始数据\太赫兹吸收测试\1mm处理后结果文件';
[num, txt, raw] = xlsread(excelPath);

% 提取时间和参考信号
t = num(:, 1);
ref_amplitude = num(:, 2);

% 提取样品信号
sampleSignals = num(:, 3:end);
sampleNames = txt(1, 3:end);

% 处理每个样品信号
for i = 1:size(sampleSignals, 2)
    sampleName = sampleNames{i};
    smp_amplitude = sampleSignals(:, i);
    
    %% 时间窗口处理
    ref_amplitude_windowed = TD_window(t, ref_amplitude, type, windowing_width);
    smp_amplitude_windowed = TD_window(t, smp_amplitude, type, windowing_width);
    
    %% 时间对齐
    [t_aligned, A_smp_pad, A_ref_pad] = time_pad(t, smp_amplitude_windowed, t, ref_amplitude_windowed);
    
    % 转置以便计算
    t_aligned = t_aligned';
    A_smp_pad = A_smp_pad';
    A_ref_pad = A_ref_pad'; 
    
    %% FFT 变换
    t_aligned = t_aligned / 1e12;  % 转换为秒
    dt = t_aligned(3) - t_aligned(2);
    fs = 1 / dt;
    f1 = (1:Zero_padding) * fs / (Zero_padding);
    f = f1' / 1e12;
    w = 2 .* pi .* f;
    
    % 参考信号 FFT
    ref_fft = fft(A_ref_pad, Zero_padding);
    ref_fft_amp = abs(ref_fft);
    ref_fft_phase = atan2(imag(ref_fft), real(ref_fft)) * 180 / pi;
    ref_fft_phase = (180 / pi) * unwrap(ref_fft_phase * pi / 180);
    
    % 样本信号 FFT
    smp_fft = fft(A_smp_pad, Zero_padding);
    smp_fft_amp = abs(smp_fft);
    smp_fft_phase = atan2(imag(smp_fft), real(smp_fft)) * 180 / pi;
    smp_fft_phase = (180 / pi) * unwrap(smp_fft_phase * pi / 180);
    
    %% 计算透射率和电磁屏蔽效能
    PowerTrans = 100 * ((smp_fft_amp.^2) ./ (ref_fft_amp.^2));
    EMISE = -10 * log10(PowerTrans / 100);
    
    %% 计算折射率
    dphase = ref_fft_phase - smp_fft_phase;
    refrac_index = 1 + dphase .* c0 ./ (180 / pi) ./ (smp_thickness .* w);
    refrac_index = refrac_index(Num_eff_0:Num_eff_1, :);
    
    %% 计算折射率和吸收系数
    dphase = abs(ref_fft_phase - smp_fft_phase);
    r_amp = smp_fft_amp ./ ref_fft_amp;
    n = 1 + (dphase .* c0 ./ (180 / pi) ./ (smp_thickness .* w));
    k = c0 ./ (smp_thickness .* w) .* log((n .* (1 + nn).^2) ./ (r_amp .* (nn + n).^2));
    absoe = 2 .* w .* k ./ c0 * 10^4;
    
                                   

    %% 计算介电常数
    diele_real = n.^2 - k.^2;
    diele_image = 2 * n .* k;
    diele_real = diele_real(Num_eff_0:Num_eff_1, :);
    diele_image = diele_image(Num_eff_0:Num_eff_1, :);
    
    %% 计算电导率
    tf_spec = smp_fft ./ ref_fft;
    cond = conj(((nn + 1) / (Z0 * smp_thickness * 1e-6)) .* ((1 ./ (tf_spec)) - 1));
    %cond_image = imag(cond(Num_eff_0:Num_eff_1, :));  % 修改这行
   % cond_real = real(cond(Num_eff_0:Num_eff_1, :));   % 修改这行
    size(diele_image)
    size(w)
    size(eps0)
    f=f(Num_eff_0:Num_eff_1, :)
    w=w(Num_eff_0:Num_eff_1, :);
    %% permittivity of bulk（块体材料适用）
    cond_real= diele_image.*w.*eps0.*1e6;           %电导率实部    [S/m]
    cond_image = (1-diele_real).*w.*eps0.*1e6;        %电导率虚部
    %Z = sqrt(1./eps_b).*tanh(1i*w.*smp_thickness./c0.*sqrt(eps_b));   %阻抗匹配系数
    con_all = [cond_real, cond_image];
    % 减少数据点数
    step = 3;
    eps_b1=diele_real
    eps_b2=diele_image
    cond_b1=cond_real
    cond_b2= cond_image
    class(eps_b1)
    f_reduced = f(1:step:end);
    eps_b1_reduced = eps_b1(1:step:end);
    eps_b2_reduced = eps_b2(1:step:end);
    cond_b1_reduced = cond_b1(1:step:end);
    cond_b2_reduced = cond_b2(1:step:end);
    size(diele_image)
    size(w)
    size(eps0)
    size(f_reduced)

    % 创建表格
    data_table = table(f_reduced, eps_b1_reduced, eps_b2_reduced, cond_b1_reduced, cond_b2_reduced, ...
        'VariableNames', {'Frequency_THz', 'Epsilon_Real', 'Epsilon_Imag', 'Conductivity_Real_S_m', 'Conductivity_Imag_S_m'});
    
    % 保存为Excel文件
    output_filename = sprintf('%s_Electromagnetic_Properties.xlsx', sampleName);
    writetable(data_table, fullfile(output_folder, output_filename));
    

    %% 提取有效部分
    PowerTrans = PowerTrans(Num_eff_0:Num_eff_1, :);
    EMISE = EMISE(Num_eff_0:Num_eff_1, :);
    con_all = con_all
    f_eff = f;
    disp('Dimensions of variables:');
    disp(['f_eff: ', num2str(size(f_eff))]);
    disp(['PowerTrans: ', num2str(size(PowerTrans))]);
    disp(['EMISE: ', num2str(size(EMISE))]);
    disp(['diele_real: ', num2str(size(diele_real))]);
    disp(['diele_image: ', num2str(size(diele_image))]);
    disp(['cond_real: ', num2str(size(cond_real))]);
    disp(['cond_image: ', num2str(size(cond_image))]);

    % 存储结果
    results = [f_eff, PowerTrans, EMISE, diele_real, diele_image, cond_real, cond_image];
    
    % 创建列名
    columnNames = {'Frequency', 'PowerTrans', 'EMISE', 'Diele_Real', 'Diele_Image', 'Cond_Real', 'Cond_Image'};
    
    % 将列名和结果合并
    resultsWithHeaders = [columnNames; num2cell(results)];
    
    % 保存结果到文件
    outputFileName = sprintf('%s-Results.csv', sampleName);
    outputFilePath = fullfile('D:\张雨林\导电力学协同优化\原始数据\太赫兹吸收测试\1mm处理后结果文件', outputFileName);
    writecell(resultsWithHeaders, outputFilePath);
end




