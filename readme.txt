ER/NR主带与展宽对比程序说明
==========================

本包用于复现我们的ER/NR主带及展宽对比图，保存了本次使用的本地NEST源码、参数、输入和统计程序，供与LZ生产流程对照。
以下文件路径均相对于解压后的代码包根目录。

一、程序逻辑

平坦能谱输入 → 调用NEST生成逐事件的S1c、S2c → 按固定S1c切片 → 计算主线和展宽 → 与论文曲线叠图。

ER、NR的S1c切片全宽分别为20、10 phd。在每个切片内，对Y=log10(S2c/phd)计算均值、10%和90%分位数，以及均值±标准差（ddof=0）。目前不额外乘经验展宽系数。10%—90%对应中央80%分位区间；这里的标准差是Y的标准差，不是复合涨落参数sigma_p。

二、建议重点审阅的文件

1. response/src/nestDiamxCore.cpp
   组织能量和位置采样、探测器配置及NEST响应调用，输出逐事件的S1c、S2c。

2. third_party/nest/include/Detectors/LZ_WS2024.hh
   实际使用的探测器、产额和涨落参数。文件名LZ_WS2024.hh和函数名SetEFT2024Configuration()沿用了2024的命名，但本流程实际调用该函数，启用的是与LZ2026论文对应的EFT校准配置，而非文件中的旧WS2024默认配置。该调用会切换相关增益设置及ER/NR产额、涨落参数向量；请以实际启用的配置和参数数值为准，不以名称中的年份判断参数版本。

3. third_party/nest/src/NEST.cpp
   实际调用的NEST核心公式，包括产额、涨落及S1/S2响应计算。

4. compute_bands.py
   保留有限且为正的S1c/S2c事件对，按固定S1c切片，计算均值、标准差和10%—90%分位数。

5. plot_bands.py
   读取统计结果和论文参考曲线，绘制主线、展宽边界及宽度比。

response_settings.json是便于审阅的参数汇总表，包含g1、g2及产额、涨落参数等。运行程序不读取此文件，单独修改它不会改变模拟。
reproduce.py是总入口，负责调度编译、事件生成、统计和绘图。

三、从头运行

解压后进入代码包根目录。需要Python及依赖；完整事件生成还需要CMake和支持C++17的编译器。

先安装Python依赖：
    python -m pip install -r requirements.txt

在已具备编译环境的Linux下运行：
    python reproduce.py --simulate

Windows使用已配置好编译环境的WSL Ubuntu-22.04时运行：
    python reproduce.py --simulate --wsl Ubuntu-22.04

上述命令依次完成：编译响应程序 → 生成ER 100万、NR 300万事件 → 计算主线和展宽 → 绘图。
新事件、统计结果、图和运行日志均保存到output/目录。

四、只重画已有结果

    python reproduce.py

此命令使用包内已有统计数据重画图，不重新生成事件，也不重新计算分位数。

