# 三体运动数值验证与八字轨道复现

三体运动数值验证与八字轨道复现

作者：____________

（单位：____________　专业：____________　学号：____________）

摘要：为检验引力轨道数值解的可信度，本文研究固定中心力、二体及三体初值问题，并复现等质量三体的八字形周期轨道。采用自编四阶龙格—库塔法，通过解析解、守恒量、步长减半、时间反演及独立 DOP853 参考轨迹交叉验证。结果表明，给定单体算例为抛物线；二体相对运动为椭圆，质心作匀速运动；三体在给定时段内出现明显分离，但不足以据此断言永久逃逸或混沌。基本算例与独立参考的最大状态差为 10⁻¹³ 至 10⁻¹² 量级。八字轨道的一周期闭合残差约为 3.90×10⁻⁸，20 周期计算仍保持轨道形态。步长扫描表明，该闭合残差主要受初值与周期有效位数限制。研究强调，应同时考察轨迹误差、守恒量及模型假设，避免仅凭图形或能量守恒判定解的正确性。

关键词：三体问题；龙格—库塔法；数值验证；八字轨道；初值敏感性

Numerical Validation of Three Body Motion and
Reproduction of the Figure Eight Orbit

Author: ____________

(Affiliation: ____________; Major: ____________; Student ID: ____________)

Abstract: This study examines the numerical accuracy of planar gravitational trajectories, progressing from motion around a fixed central mass to interacting two-body and three-body systems. A classical fourth-order Runge–Kutta method is implemented in Julia. Its results are assessed using an analytical parabolic solution, conservation laws, step refinement, time reversal, and independently implemented DOP853 reference trajectories. For the prescribed initial conditions, the fixed-center example follows a parabola, while the relative two-body motion is elliptical and its center of mass drifts uniformly. The unequal-mass three-body example exhibits substantial separation during the simulated interval, but this observation alone does not establish permanent escape or chaotic dynamics. Maximum state differences from the independent reference are approximately 10⁻¹³ to 10⁻¹² for the basic examples. Velocity perturbations are examined at two amplitudes after removing center-of-mass motion, separating physical sensitivity from numerical discretization error. The additional experiment reproduces the equal-mass figure-eight solution originally discovered numerically by Moore. Periodicity, cyclic exchange after one third of a period, minimum pairwise separation, and conservation errors are evaluated over one and twenty periods. Refining the time step reduces the difference from the reference trajectory, whereas the one-period closure residual approaches approximately 3.90×10⁻⁸. This plateau is consistent with the finite precision of the published initial state and approximate period. The results demonstrate a reproducible verification procedure while distinguishing numerical agreement from mathematical proof of periodicity, stability, or chaos. Complete source files, reference data, and validation tests accompany the study.

Keywords: three-body problem; Runge–Kutta method; numerical validation; figure-eight orbit; initial-condition sensitivity

引言

引力多体运动同时包含可解析检验的简单情形与难以获得通用显式解的复杂情形，适合用于研究数值积分误差。Moore 于 1993 年通过作用量极小化的数值探索发现了包括八字轨道在内的平面多体运动[1]；Chenciner 与 Montgomery 后来给出了等质量三体八字解的存在性证明[2]。后续稳定性与分岔研究表明，周期轨道的动力学性质需要专门分析[3]。因此，画出闭合曲线并不能替代正确性或稳定性论证。

已有工作还提供了结构保持积分方法[4]、开源多体程序 REBOUND[5]、高阶积分器 IAS15[6]及 Julia 微分方程生态[7]。这些资源为后续交叉验证提供条件。本文以课堂初值为对象，先完成单体、二体、三体基本任务，再开展八字解复现；REBOUND、IAS15 与 DifferentialEquations.jl 仅作为相关工作，本次独立对照实际使用 SciPy DOP853。

1  模型与数值方法

1.1  无量纲模型与初值

所有算例采用无量纲量，取引力常数 G=1，不引入引力软化。单体问题指中心天体固定、卫星反作用可忽略的模型，并补充取中心质量 M=1。对相互作用的 N 个质点，令状态 z 依次排列所有位置与所有速度，运动方程为



ṙᵢ = vᵢ，　v̇ᵢ = G ∑ⱼ≠ᵢ mⱼ (rⱼ − rᵢ) / ‖rⱼ − rᵢ‖³ (1-1)

两体距离趋于零时原方程奇异。本次给定时间区间内均未发生碰撞，结论不延伸到未计算时段。表 1-1 列出基本算例；位置与速度均按物体编号列出。

表1-1  基本算例的无量纲初值与时段



表1-1  基本算例的无量纲初值与时段
算例 | 质量 | 初始位置 | 初始速度 | 终止时间
单体 | 中心 M=1 | (0, 2) | (1, 0) | 50
二体 | 0.3；0.03 | (2, 2)；
(0, 0) | (0.2, −0.2)；
(−0.01, 0.01) | 20
三体 | 0.3；0.03；0.03 | (2, 2)；(0, 0)；
(−2, −2) | (0.2, −0.2)；
(0, 0)；(−0.2, 0.2) | 30

1.2  积分算法与伪代码

主程序采用 Julia[8] 编写，以经典 RK4 为基本积分器，基本任务步长 h=0.001。独立 Python 实现采用 NumPy 数组[9]与 SciPy[10]，通过 DOP853 的高阶自适应积分生成参考数据[11]；该算法的原始代码可由 Hairer 的软件页面追溯[12]。项目同时保留 Dormand–Prince 自适应方法作为扩展对照[13]，误差控制设计参照嵌入式 Runge–Kutta 方法的相关讨论[14]。

伪代码1  RK4 轨道积分

输入：质量数组 m，初始状态 z₀，终止时间 T，最大步长 h。

初始化：t ← 0，z ← z₀；保存初始状态。

当 t < T 时：

　Δt ← min(h, T − t)。

　k₁ ← f(z, t)；k₂ ← f(z + Δt k₁/2, t + Δt/2)。

　k₃ ← f(z + Δt k₂/2, t + Δt/2)。

　k₄ ← f(z + Δt k₃, t + Δt)。

　z ← z + Δt (k₁ + 2k₂ + 2k₃ + k₄)/6。

　t ← t + Δt；保存状态并检查数值是否有限。

输出：时间序列、位置速度序列及诊断量。

力函数按物体对计算距离与引力项，并组装式（1-1）。RK4 在光滑且远离碰撞的情形下具有四阶全局收敛性；当步长已足够小时，继续减小步长可能进入舍入误差主导区，不能要求误差始终按 16 倍下降。

1.3  正确性检验指标

对二体、三体计算能量 E、总动量 P、角动量 Lz 和质心 R。令 M 为总质量，rij 为物体间距离，内部能量 Ec 去除整体平动动能：



E = ½∑ᵢ mᵢ‖vᵢ‖² − ∑ᵢ<ⱼ Gmᵢmⱼ/rᵢⱼ，　Eᵢₙₜ = E − ‖P‖²/(2M) (1-2)



P = ∑ᵢ mᵢvᵢ，　Lz = ∑ᵢ mᵢ(xᵢvᵧᵢ − yᵢvₓᵢ) (1-3)



R(t) = R(0) + P(0)t/M (1-4)

参考差定义为相同采样时刻、所有无量纲状态分量绝对差的最大值。基本任务参考采样间隔为 0.1；守恒量在主积分的全部步上检查。DOP853 使用 rtol=3×10⁻¹⁴、atol=3×10⁻¹⁶，并与较宽容差结果比较。参考解仍有离散与舍入误差，因此两种积分器的一致性不构成严格误差上界。

伪代码2  多重数值验证

输入：算例集合 C；步长 h；公共采样网格 Γ；扰动集合 {10⁻⁶, 10⁻⁸}
输出：验证记录 Q（参考差、收敛差、守恒量、反演和轨道诊断）

01  Q ← ∅
02  for c ∈ C do
03      Z_h ← RK4(c, h)
04      Z_h/2 ← RK4(c, h/2)
05      Q[c].step ← max_{t∈Γ} ‖Z_h(t) − Z_h/2(t)‖∞
06      if c 有解析解 A_c(t) then
07          Q[c].reference ← max_{t∈Γ} ‖Z_h(t) − A_c(t)‖∞
08      else
09          Z_ref ← DOP853(c, Γ)
10          Q[c].reference ← max_{t∈Γ} ‖Z_h(t) − Z_ref(t)‖∞
11      end if
12      if c 为相互作用的多体系统 then
13          Q[c].conservation ← 漂移(E, P, Lz, R(t) − R(0) − P(0)t/M)
14      end if
15      c_rev ← 将 c 的初态替换为 (r_h(T), −v_h(T))
16      Z_rev ← RK4(c_rev, h)
17      Q[c].reverse ← ‖(r_rev(T), −v_rev(T)) − z₀‖∞
18      if c 配置了初值扰动 then
19          for δ ∈ {10⁻⁶, 10⁻⁸} do
20              Z_δ,h ← RK4(扰动(c, δ), h)
21              Z_δ,h/2 ← RK4(扰动(c, δ), h/2)
22              D_δ,h(t) ← 去质心位置差(Z_δ,h(t), Z_h(t))
23              D_δ,h/2(t) ← 去质心位置差(Z_δ,h/2(t), Z_h/2(t))
24              Q[c,δ].response ← D_δ,h(T)
25              Q[c,δ].step ← max_{t∈Γ} |D_δ,h(t) − D_δ,h/2(t)|
26          end for
27      end if
28      if c 为八字算例 then
29          Q[c].closure ← ‖Z_h(T) − Z_h(0)‖∞
30          Q[c].symmetry ← max_t ‖Z_h(t + T/3) − ΠZ_h(t)‖∞
31          Q[c].distance ← min_{t, i<j} ‖rᵢ(t) − rⱼ(t)‖₂
32      end if
33  end for
34  return Q

其中，去质心位置差是两条轨迹在各自质心系中所有物体位置差的欧氏范数；Π 表示三体编号的循环置换。

2  基本任务结果与分析

2.1  固定中心力的抛物线解

由初值可得单位质量机械能为 v²/2−1/r=0，角动量为 −2，故轨道为抛物线而非椭圆。设辅助变量 D 满足 D+D³/3=t/4，则解析位置与速度为



D = 2sinh[arsinh(3t/8)/3]，　x = 4D，　y = 2 − 2D² (2-1)



vₓ = 1/(1+D²)，　vᵧ = −D/(1+D²) (2-2)

消去 D 得 x²=16−8y。RK4 与解析完整状态的最大误差约为 4.14×10⁻¹³，轨道形态和误差检验一致。该结论依赖补充假设 GM=1；若中心质量不同，同一初始速度可能对应不同类型的圆锥曲线。



图2-1  固定中心力算例的抛物线轨道

[图] 图2-1  固定中心力算例的抛物线轨道

2.2  二体相对轨道与质心漂移

本算例初始总动量为 (0.0597, −0.0597)，故惯性系轨道包含质心漂移。去除质心平动后，相对运动的半长轴约为 2.273585、偏心率约为 0.244039、周期约为 37.4964。计算终点 t=20 尚未覆盖一整周，不能将未闭合的绘图误认为算法失败。虽然总能量约为 0.00882102，但内部能量为 −0.00197925，二体相对运动仍受束缚。



图2-2  二体运动的质心系轨迹

[图] 图2-2  二体运动的质心系轨迹

2.3  三体运动的有限时段结论

三体初始质心为 (1.5, 1.5)，质心速度为 (0.15, −0.15)。总能量约为 0.00810883，内部能量约为 8.83118×10⁻⁶。在 t=30 时，三组两体距离约为 2.5783、14.5385 和 15.3203，表明第三体在观测时段内远离前两体。全步最小距离约为 1.4621。有限时段的分离现象尚不能证明永久逃逸，也不能直接证明混沌。



图2-3  三体运动的质心系轨迹

[图] 图2-3  三体运动的质心系轨迹

2.4  误差收敛与微扰响应

表2-1  基本算例的完整状态误差



表2-1  基本算例的完整状态误差
算例 | 与 DOP853 最大差 | 步长减半最大差 | 时间反演误差
单体 | 1.51×10⁻¹³ | 2.31×10⁻¹³ | 1.33×10⁻¹⁵
二体 | 1.57×10⁻¹³ | 1.98×10⁻¹³ | 6.76×10⁻¹⁵
三体 | 1.47×10⁻¹² | 1.34×10⁻¹³ | 5.31×10⁻¹⁵

表 2-1 使用 h=0.001，减半对照为 h=0.0005。较粗步长扫描呈现 RK4 的收敛趋势，细步长结果达到参考与浮点误差平台。二体全步能量、动量、角动量、质心误差分别约为 9.54×10⁻¹⁷、4.30×10⁻¹⁶、1.39×10⁻¹⁵ 和 3.20×10⁻¹⁴；三体对应约为 3.59×10⁻¹⁶、8.19×10⁻¹⁶、3.44×10⁻¹⁵ 和 3.91×10⁻¹⁴。守恒检验支持数值一致性，但必须结合轨迹对照使用。



图2-4  基本算例对独立参考轨迹的步长收敛

[图] 图2-4  基本算例对独立参考轨迹的步长收敛

二体对第二体初始 y 方向速度加扰动，三体对第三体初始 x 方向速度加扰动，分别取 10⁻⁶ 和 10⁻⁸。比较量为去除各自质心后的全体位置差的欧氏范数，不能与单个物体位置差混用。

表2-2  去质心后的末态位置微扰响应



表2-2  去质心后的末态位置微扰响应
算例 | 速度扰动 10⁻⁶ | 速度扰动 10⁻⁸ | 终止时间
二体 | 2.619520×10⁻⁵ | 2.619521×10⁻⁷ | 20
三体 | 3.398810×10⁻⁵ | 3.398811×10⁻⁷ | 30

两档扰动下响应近似随扰动幅值线性缩放。扰动信号的步长减半差处于约 10⁻¹³ 量级，明显低于表 2-2 的响应。因而本实验可以分辨有限时间内的初值敏感性；没有计算收敛的最大 Lyapunov 指数，不能把这些响应称为已证实的指数发散。

3  八字周期轨道复现

3.1  文献对应与初始条件

本节复现 Moore 所发现的 figure-eight solution，即三个等质量质点沿同一八字曲线依次运动的解[1-2]，不复现其寻找轨道的作用量极小化算法。数值初值采用 GeometricProblems.jl 的公开三体示例[15]并按本文编号重排。该示例是具体数值的来源，不能把这些小数直接归为 Moore 论文给出的参数。

取 G=1、m₁=m₂=m₃=1、近似周期 T=6.32591398。它是与基本任务质量比例不同的独立算例。表 3-1 的初态满足质心在原点、总动量和总角动量为零。

表3-1  八字轨道的无量纲初始条件



表3-1  八字轨道的无量纲初始条件
物体 | 初始位置 (x, y) | 初始速度 (vx, vy)
1 | (−0.97000436, 0.24308753) | (0.466203685, 0.432365730)
2 | (0, 0) | (−0.932407370, −0.864731460)
3 | (0.97000436, −0.24308753) | (0.466203685, 0.432365730)



图3-1  八字轨道及 t=T/12 时的三体位置

[图] 图3-1  八字轨道及 t=T/12 时的三体位置

3.2  周期闭合与舞蹈对称性

按每周期 300 至 9600 步扫描。周期闭合残差为末态与初态的最大分量差，参考差为相同初值的 RK4 与 DOP853 在公共网格上的最大分量差，两者检验对象不同。

表3-2  一周期闭合与独立参考差



表3-2  一周期闭合与独立参考差
每周期步数 | 完整状态闭合残差 | 与 DOP853 最大差
300 | 1.58×10⁻⁶ | 2.61×10⁻⁶
600 | 4.31×10⁻⁸ | 1.44×10⁻⁷
1200 | 3.55×10⁻⁸ | 8.37×10⁻⁹
2400 | 3.88×10⁻⁸ | 5.04×10⁻¹⁰
4800 | 3.90×10⁻⁸ | 3.05×10⁻¹¹
9600 | 3.90×10⁻⁸ | 1.62×10⁻¹²

细化步长后参考差持续降低，而闭合残差停留在约 3.90×10⁻⁸；独立参考轨迹自身的一周期闭合残差也约为 3.90×10⁻⁸。因此该平台与给定初值、周期的小数截断一致，不能解释为 RK4 不收敛。每周期 600 步处闭合残差偶然较小也不代表更准确。

在本文物体顺序下，时间推进 T/3 后，完整状态应对应原状态的物体置换 (1, 2, 3)→(3, 1, 2)。采用每周期 2400 步，第一周期内的位置与速度联合置换残差约为 2.83×10⁻⁸，支持三体沿同一路径依次追逐的舞蹈对称性。

3.3  二十周期验证与局限

取 h=T/2400 积分 20 周期，位置闭合残差约为 5.75×10⁻⁷，速度闭合残差约为 7.49×10⁻⁷；相同初值下与独立参考的最大状态差约为 2.24×10⁻⁸。全程能量绝对漂移不超过约 8.62×10⁻¹¹，最小两体距离约为 0.69053。收紧参考容差引起的轨迹差约为 5.22×10⁻⁹，说明长时间对照亦应报告参考自身的不确定性。



图3-2  八字轨道的步长扫描与二十周期检验

[图] 图3-2  八字轨道的步长扫描与二十周期检验

图 3-2 同时展示积分误差、周期闭合和局部相位投影。局部相位估计将位置误差投影到参考速度方向，仅用于区分沿轨迹与横向偏离，并非全局相位的精确重建。20 周期内保持轨道形态不等于证明无限时间稳定性；严格的轨道存在性与稳定性论证应与文献中的变分理论和特征乘子分析区分[2-3]。

4  可复现性与项目组织

源代码、独立参考数据、图像和检查脚本汇集于独立 GitHub 仓库：

https://github.com/zhangzhiyang2006/three-body-numerical-study

仓库提供不含运行输出的 Notebook、可维护源文件、Julia 环境声明、Python 依赖文件及参考数据的 SHA-256 校验值。主程序生成图表，参考脚本以独立的成对引力实现生成 DOP853 数据。数值结果保存为 TOML 文件，便于逐项核对。本次独立项目通过 23 项 Julia 数值检查及 3 项 Python 交付一致性检查。具体安装与执行步骤见仓库 README。

源文件保留相对路径；姓名、学号、个人邮箱、本机绝对路径、临时日志与其他课程作业不纳入发布内容。作者信息在本地提交稿中自行填写，GitHub 仓库只保留匿名版本。仓库访问受其可见性设置控制。

5  总结与展望

本文通过由简到繁的算例构建了可重复的轨道验证流程。固定中心力算例与抛物线解析解吻合；二体相对轨道及质心运动满足理论预期；三体计算在给定时段内与独立积分结果一致，并可分辨两档微小速度扰动的响应。在此基础上，等质量八字轨道的周期闭合、T/3 状态置换和 20 周期轨迹对照均支持复现结果。

本研究仍采用有限精度初值和有限积分时间。后续可用高精度射击法联合修正初值与周期，降低闭合残差平台；通过变分方程与特征乘子研究局部稳定性；在相同误差预算下比较 RK4、辛积分与 IAS15 的计算成本，并为近碰撞情形引入事件检测和适当正则化。上述工作属于后续研究，本文不将其作为已完成结果。

参考文献

[1] MOORE C. Braids in classical dynamics[J]. Physical Review Letters, 1993, 70(24): 3675-3679. DOI: 10.1103/PhysRevLett.70.3675.

[2] CHENCINER A, MONTGOMERY R. A remarkable periodic solution of the three-body problem in the case of equal masses[J]. Annals of Mathematics, 2000, 152(3): 881-901. DOI: 10.2307/2661357.

[3] GALÁN J, MUÑOZ-ALMARAZ F J, FREIRE E, et al. Stability and bifurcations of the figure-8 solution of the three-body problem[J]. Physical Review Letters, 2002, 88(24): 241101. DOI: 10.1103/PhysRevLett.88.241101.

[4] HAIRER E, LUBICH C, WANNER G. Geometric numerical integration illustrated by the Störmer–Verlet method[J]. Acta Numerica, 2003, 12: 399-450. DOI: 10.1017/S0962492902000144.

[5] REIN H, LIU S F. REBOUND: An open-source multi-purpose N-body code for collisional dynamics[J]. Astronomy & Astrophysics, 2012, 537: A128. DOI: 10.1051/0004-6361/201118085.

[6] REIN H, SPIEGEL D S. IAS15: A fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine precision over a billion orbits[J]. Monthly Notices of the Royal Astronomical Society, 2015, 446(2): 1424-1437. DOI: 10.1093/mnras/stu2164.

[7] RACKAUCKAS C, NIE Q. DifferentialEquations.jl—A performant and feature-rich ecosystem for solving differential equations in Julia[J]. Journal of Open Research Software, 2017, 5(1): 15. DOI: 10.5334/jors.151.

[8] BEZANSON J, EDELMAN A, KARPINSKI S, et al. Julia: A fresh approach to numerical computing[J]. SIAM Review, 2017, 59(1): 65-98. DOI: 10.1137/141000671.

[9] HARRIS C R, MILLMAN K J, VAN DER WALT S J, et al. Array programming with NumPy[J]. Nature, 2020, 585: 357-362. DOI: 10.1038/s41586-020-2649-2.

[10] VIRTANEN P, GOMMERS R, OLIPHANT T E, et al. SciPy 1.0: Fundamental algorithms for scientific computing in Python[J]. Nature Methods, 2020, 17: 261-272. DOI: 10.1038/s41592-019-0686-2.

[11] SCIPY DEVELOPERS. scipy.integrate.DOP853[EB/OL]. [2026-10-01]. https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.DOP853.html.

[12] HAIRER E. Fortran codes for differential equations[EB/OL]. [2026-10-01]. https://www.unige.ch/~hairer/software.html.

[13] DORMAND J R, PRINCE P J. A family of embedded Runge-Kutta formulae[J]. Journal of Computational and Applied Mathematics, 1980, 6(1): 19-26. DOI: 10.1016/0771-050X(80)90013-3.

[14] SHAMPINE L W. Some practical Runge-Kutta formulas[J]. Mathematics of Computation, 1986, 46(173): 135-150. DOI: 10.1090/S0025-5718-1986-0815836-3.

[15] JULIAGNI CONTRIBUTORS. Three body problem: GeometricProblems.jl[EB/OL]. [2026-10-01]. https://juliagni.github.io/GeometricProblems.jl/latest/three_body_problem/.
