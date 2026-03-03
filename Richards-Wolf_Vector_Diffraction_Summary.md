# Richards-Wolf 矢量衍射理论：原理、推导与应用详解

**目录**

1.  [引言：从标量到矢量](#1-引言从标量到矢量)
2.  [理论基础与物理模型](#2-理论基础与物理模型)
    *   [2.1 几何坐标系定义](#21-几何坐标系定义)
    *   [2.2 Debye-Wolf 积分原理](#22-debye-wolf-积分原理)
    *   [2.3 Debye 近似的物理假设与菲涅尔数验证](#23-debye-近似的物理假设与菲涅尔数验证)
3.  [Richards-Wolf 积分公式的完整推导](#3-richards-wolf-积分公式的完整推导)
    *   [3.1 亥姆霍兹方程与角谱表示](#31-亥姆霍兹方程与角谱表示)
    *   [3.2 步骤一：切趾因子 (Apodization) 的推导](#32-步骤一切趾因子-apodization-的推导)
    *   [3.3 步骤二：偏振变换矩阵 $\mathcal{M}$ 的推导](#33-步骤二偏振变换矩阵-mathcalm-的推导)
    *   [3.4 步骤三：矢量强度因子 $\mathbf{a}$ 的构建](#34-步骤三矢量强度因子-mathbfa-的构建)
    *   [3.5 最终积分形式](#35-最终积分形式)
4.  [解析解求解：贝塞尔函数展开](#4-解析解求解贝塞尔函数展开)
    *   [4.1 Jacobi-Anger 展开技术](#41-jacobi-anger-展开技术)
    *   [4.2 线偏振光的场分量推导 ($I_0, I_1, I_2$)](#42-线偏振光的场分量推导-i_0-i_1-i_2)
    *   [4.3 任意偏振态的通用计算](#43-任意偏振态的通用计算)
5.  [物理意义深度解析](#5-物理意义深度解析)
    *   [5.1 矢量效应与标量理论的对比](#51-矢量效应与标量理论的对比)
    *   [5.2 纵向分量 $E_z$ 与去偏振效应](#52-纵向分量-e_z-与去偏振效应)
    *   [5.3 相位异常与古依相移 (Gouy Phase Shift)](#53-相位异常与古依相移-gouy-phase-shift)
6.  [像差建模与广义光瞳函数](#6-像差建模与广义光瞳函数)
    *   [6.1 广义光瞳函数](#61-广义光瞳函数)
    *   [6.2 数值计算策略：FFT 方法](#62-数值计算策略fft-方法)
7.  [总结与应用领域](#7-总结与应用领域)

---

## 1. 引言：从标量到矢量

在光学系统中，当数值孔径（Numerical Aperture, NA）较低时，光线的会聚角较小，电磁场的标量衍射理论（如 Fresnel 或 Fraunhofer 衍射）足以精确描述焦斑分布。此时，我们将光视为标量场 $U$，忽略电场矢量的方向性。

然而，当 **$NA > 0.7$** 时，光线会以极大的角度（半角 $\alpha > 45^\circ$）会聚。此时产生两个不可忽略的效应：
1.  **矢量特性 (Vectorial Nature)**：光是一种横波，不同方向会聚的光线，其电场矢量方向不同，简单的代数叠加失效，必须进行矢量叠加。
2.  **偏振效应 (Polarization Effects)**：原本沿横向（x/y）振动的光，经大角度折射后，会产生显著的纵向分量（$z$ 分量）。

为了解决这一问题，Peter Debye 提出了基于平面波谱叠加的积分表示，后经 E. Wolf 和 B. Richards 完善，形成了描述无像差高 NA 聚焦场的核心理论——**Richards-Wolf 矢量衍射理论**。它不是对标量理论的修正，而是直接基于麦克斯韦方程组的严谨解。

---

## 2. 理论基础与物理模型

### 2.1 几何坐标系定义

我们构建如下物理模型：
*   **光学系统**：一个高数值孔径的无像差（Aplanatic）物镜，焦距为 $f$，最大孔径半角为 $\alpha$。
    $$ \text{NA} = n \sin \alpha $$
*   **坐标系**：
    *   原点 $O$ 位于几何焦点。
    *   **球坐标 $(r, \theta, \phi)$**: 描述入瞳及参考球面上的波前。$\theta$ 为极角（相对于光轴 $z$），$\phi$ 为方位角。
    *   **柱坐标 $(\rho, \varphi, z)$**: 描述像方焦平面附近的观察点 $P$ 位置。
        $$ x_p = \rho \cos\varphi, \quad y_p = \rho \sin\varphi, \quad z_p = z $$
*   **参考球面 (Reference Sphere)**：以焦点 $O$ 为球心，半径为 $f$ 的球面，代表出瞳处的理想波前。

### 2.2 Debye-Wolf 积分原理
**Debye-Wolf 积分**是一个广泛的数学框架，它指出：**焦点附近的场可以看作是来自出瞳孔径的无数个平面波的相干叠加**。它本质上是一种角谱（Angular Spectrum）表示法。

**Richards-Wolf 理论**则是 Debye-Wolf 积分在“无像差高数值孔径系统”中的具体实现。它特指引入了以下三个物理限制后的积分形式：
1.  **阿贝正弦条件 (Abbe Sine Condition)**。
2.  **矢量偏振传输 (Vector Polarization)**。
3.  **能量守恒切趾 (Apodization)**。

### 2.3 Debye 近似的物理假设与菲涅尔数验证
Debye 积分不仅仅是基尔霍夫衍射公式的简化，它基于以下核心物理假设：
1.  **平面波叠加**：在焦点附近的小区域内，弯曲的波前可以分解为一个个局部的平面波分量。
2.  **孔径限制**：积分区域严格限制在孔径角 $\Omega$ 内。
3.  **忽略边缘衍射 (Edge Diffraction)**：假设孔径光阑边缘产生的衍射波对焦点的贡献可忽略。

**为什么可以忽略边缘衍射？（菲涅尔数判据）**
菲涅尔数 $N$ 描述了孔径边缘衍射波与直接几何波在焦点处的相位差（以半波长为单位）：
$$ N = \frac{a^2}{\lambda f} $$
其中 $a$ 是孔径半径。

**典型计算案例**：
假设使用高倍油浸物镜：100x, NA=1.4, 工作波长 $\lambda=0.5 \mu m$, 焦距 $f \approx 1.6 mm$。
$$ a \approx f \cdot \text{NA} \approx 2.24 mm $$
$$ N = \frac{(2.24 \times 10^{-3} m)^2}{0.5 \times 10^{-6} m \times 1.6 \times 10^{-3} m} \approx 6272 $$

*   当 $N$ 较小（如针孔相机，$N<10$）时，边缘衍射会导致显著的**焦移 (Focal Shift)**，即实际最大光强点不在几何焦点。
*   当 $N \gg 1$ （如显微物镜，$N \sim 1000$）时，边缘效应微乎其微，几何焦点与波动光学焦点重合，Debye 近似极其精确。

---

## 3. Richards-Wolf 积分公式的完整推导

### 3.1 亥姆霍兹方程与角谱表示
根据麦克斯韦方程组，像方空间的电场 $\mathbf{E}(P)$ 可以表示为加权平面波的叠加：

$$ \mathbf{E}(P) = -\frac{ikf}{2\pi} \iint_{\Omega} \mathbf{A}(\theta, \phi) e^{i \mathbf{k} \cdot \mathbf{r}_P} d\Omega $$

其中：
*   **$\mathbf{k}$ (波矢量)**：方向由 $(\theta, \phi)$ 决定。
    $$ \mathbf{k} = k (\sin\theta \cos\phi, \sin\theta \sin\phi, \cos\theta) $$
*   **$\mathbf{r}_P$ (位置矢量)**：观察点坐标。
*   **$\mathbf{A}(\theta, \phi)$**: 矢量角谱，包含了振幅和偏振信息。
*   **$d\Omega$**: 立体角元 $\sin\theta d\theta d\phi$。

相位的展开形式为：
$$ \mathbf{k} \cdot \mathbf{r}_P = k (x \sin\theta \cos\phi + y \sin\theta \sin\phi + z \cos\theta) $$
$$ = k z \cos\theta + k \rho \sin\theta \cos(\phi - \varphi) $$

接下来的核心任务是构建 **$\mathbf{A}(\theta, \phi)$**，这需要两步修正。

### 3.2 步骤一：切趾因子 (Apodization) 的推导
光线经过透镜系统后，能量密度会发生改变。**切趾因子 $P(\theta)$** 描述了这种振幅变化，它来源于**能量守恒**。

**推导过程**：
1.  **入瞳能量**：考虑入瞳平面上的环形区域，半径 $r$，宽度 $dr$。假设入射光强度均匀为 $I_{in}$。
    能量通量：$d\Phi_{in} = I_{in} \cdot 2\pi r dr$
2.  **出瞳能量**：对应于参考球面上的环带，角度 $\theta$，角宽度 $d\theta$。球半径为 $f$。
    能量通量：$d\Phi_{out} = I_{out}(\theta) \cdot (\text{环带面积}) = I_{out}(\theta) \cdot (2\pi f \sin\theta) (f d\theta)$
3.  **引入阿贝正弦条件**：对于无像差物镜，满足 $r = f \sin\theta$。
    微分得：$dr = f \cos\theta d\theta$。
4.  **能量守恒联立**：令 $d\Phi_{in} = d\Phi_{out}$。
    $$ I_{in} \cdot 2\pi (f \sin\theta) (f \cos\theta d\theta) = I_{out}(\theta) \cdot 2\pi f^2 \sin\theta d\theta $$
    $$ I_{in} \cos\theta = I_{out}(\theta) $$
5.  **振幅关系**：电场振幅 $E \propto \sqrt{I}$。
    $$ P(\theta) = \frac{E_{out}}{E_{in}} = \sqrt{\frac{I_{out}}{I_{in}}} = \sqrt{\cos\theta} $$

**物理意义**：在高 NA 下，光束截面面积随角度变化，为了守恒能量，大角度光线的振幅必须按照 $\sqrt{\cos\theta}$ 衰减。这相当于一个天然的低通滤波器。

### 3.3 步骤二：偏振变换矩阵 $\mathcal{M}$ 的推导
光线经过透镜折射后，传播方向改变 $\theta$ 角。根据电磁波的横波性质，电场矢量 $\mathbf{E}$ 必须始终垂直于传播方向 $\mathbf{k}$。因此，电场矢量必须随光线一起旋转。

我们通过**坐标旋转法**推导该变换矩阵 $\mathcal{M}$：
$$ \mathbf{E}_{ref} \propto \mathcal{R}_z(\phi) \cdot \mathcal{R}_y(\theta) \cdot \mathcal{R}_z(-\phi) \cdot \mathbf{E}_{inc} $$

1.  **分解与对齐 $\mathcal{R}_z(-\phi)$**：
    将入射电场分解为**径向分量 (p)** 和 **方位角分量 (s)**。这相当于将实验室坐标系顺时针旋转 $\phi$，使当前光线位于子午面（$x-z$ 平面）。
    $$ \mathcal{R}_z(-\phi) = \begin{pmatrix} \cos\phi & \sin\phi & 0 \\ -\sin\phi & \cos\phi & 0 \\ 0 & 0 & 1 \end{pmatrix} $$
2.  **折射旋转 $\mathcal{R}_y(\theta)$**：
    在子午面内，光线向光轴弯曲 $\theta$ 角。
    *   垂直于入射面的分量 (s-pol) 方向不变。
    *   平行于入射面的分量 (p-pol) 随光线旋转 $\theta$ 角。
    $$ \mathcal{R}_y(\theta) = \begin{pmatrix} \cos\theta & 0 & \sin\theta \\ 0 & 1 & 0 \\ -\sin\theta & 0 & \cos\theta \end{pmatrix} $$
3.  **复原 $\mathcal{R}_z(\phi)$**：
    将坐标系逆时针转回原来的方位角 $\phi$。

**最终矩阵 $\mathcal{M}$**（三个矩阵相乘的结果）：
$$ \mathcal{M}(\theta, \phi) = \begin{pmatrix}
1 + (\cos\theta - 1)\cos^2\phi & (\cos\theta - 1)\cos\phi \sin\phi & -\sin\theta \cos\phi \\
(\cos\theta - 1)\cos\phi \sin\phi & 1 + (\cos\theta - 1)\sin^2\phi & -\sin\theta \sin\phi \\
\sin\theta \cos\phi & \sin\theta \sin\phi & \cos\theta
\end{pmatrix} $$

### 3.4 步骤三：矢量强度因子 $\mathbf{a}$ 的构建
综合切趾和偏振变换，像方空间的矢量因子 $\mathbf{a}(\theta, \phi)$ 为：
$$ \mathbf{a}(\theta, \phi) = P(\theta) \cdot \mathcal{M}(\theta, \phi) \cdot \mathbf{E}_{inc} $$

**案例：x-轴线偏振入射**
设入射光 $\mathbf{E}_{inc} = [1, 0, 0]^T$。 $\mathbf{a}$ 等于矩阵 $\mathcal{M}$ 的第一列乘以 $\sqrt{\cos\theta}$：
$$ \mathbf{a} = \sqrt{\cos\theta} \begin{pmatrix}
1 + (\cos\theta - 1)\cos^2\phi \\
(\cos\theta - 1)\cos\phi \sin\phi \\
\sin\theta \cos\phi
\end{pmatrix} $$
利用三角恒等式降幂处理，可将其改写为便于积分的形式：
$$ a_x = \sqrt{\cos\theta} \left[ \frac{1+\cos\theta}{2} + \frac{\cos\theta-1}{2}\cos 2\phi \right] $$
$$ a_y = \sqrt{\cos\theta} \left[ \frac{\cos\theta-1}{2}\sin 2\phi \right] $$
$$ a_z = \sqrt{\cos\theta} \left[ -\sin\theta\cos\phi \right] $$

### 3.5 最终积分形式
将 $\mathbf{a}$ 和相位因子代入原始积分：
$$ \mathbf{E}(\rho, \varphi, z) = -\frac{ikf}{2\pi} \int_{0}^{\alpha} \int_{0}^{2\pi} \mathbf{a}(\theta, \phi) e^{ikz\cos\theta} e^{ik\rho\sin\theta\cos(\phi-\varphi)} \sin\theta \, d\phi \, d\theta $$

---

## 4. 解析解求解：贝塞尔函数展开

直接计算双重积分较为困难，我们利用 **Jacobi-Anger 展开**将 $\phi$ 维度的积分转化为贝塞尔函数，从而得到半解析解。

### 4.1 Jacobi-Anger 展开技术
核心恒等式：
$$ \int_{0}^{2\pi} e^{ix\cos\beta} e^{in\beta} d\beta = 2\pi i^n J_n(x) $$
特例：
*   常数项积分 ($n=0$) $\rightarrow 2\pi J_0(x)$
*   $\cos\phi$ 项积分 ($n=1$) $\rightarrow 2\pi i J_1(x)$
*   $\cos 2\phi$ 项积分 ($n=2$) $\rightarrow -2\pi J_2(x)$

### 4.2 线偏振光的场分量推导 ($I_0, I_1, I_2$)
将 3.4 节中的 $a_x, a_y, a_z$ 代入积分，并利用上述恒等式完成 $\phi$ 积分，我们定义三个关键的径向积分 $I_0, I_1, I_2$：

$$
\begin{aligned}
I_0(\rho, z) &= \int_{0}^{\alpha} \sqrt{\cos\theta} \sin\theta (1+\cos\theta) J_0(k\rho\sin\theta) e^{ikz\cos\theta} d\theta \\
I_1(\rho, z) &= \int_{0}^{\alpha} \sqrt{\cos\theta} \sin^2\theta J_1(k\rho\sin\theta) e^{ikz\cos\theta} d\theta \\
I_2(\rho, z) &= \int_{0}^{\alpha} \sqrt{\cos\theta} \sin\theta (1-\cos\theta) J_2(k\rho\sin\theta) e^{ikz\cos\theta} d\theta
\end{aligned}
$$

**x-线偏振入射下的最终场分布**：
$$ \mathbf{E}(\rho, \varphi, z) = -\frac{ikf}{2} \begin{pmatrix}
I_0 + I_2 \cos(2\varphi) \\
I_2 \sin(2\varphi) \\
-2i I_1 \cos(\varphi)
\end{pmatrix} $$

### 4.3 任意偏振态的通用计算
Richards-Wolf 理论是一个线性理论，对于任意偏振态（如圆偏振、径向偏振），可以遵循以下通用流程：

1.  **表达输入场**: 写出入瞳处的电场矢量 $\mathbf{E}_{inc}(\phi)$。
    *   *圆偏振*: $\mathbf{E}_{inc} = \frac{1}{\sqrt{2}}[1, \pm i, 0]^T$
    *   *径向偏振*: $\mathbf{E}_{inc} = [\cos\phi, \sin\phi, 0]^T$
2.  **矩阵预处理**: 计算 $\mathbf{a}(\theta, \phi) = \sqrt{\cos\theta} \cdot \mathcal{M} \cdot \mathbf{E}_{inc}$。
    *   例如，径向偏振经矩阵乘法后，$\mathbf{a} = \sqrt{\cos\theta} [\cos\theta\cos\phi, \cos\theta\sin\phi, \sin\theta]^T$。
3.  **积分求解**:
    *   **半解析法**: 如果 $\mathbf{a}$ 能化简为 Fourier 级数形式，则继续使用贝塞尔函数求解。
    *   **数值积分法**: 如果 $\mathbf{E}_{inc}$ 很复杂，直接对 $\theta, \phi$ 网格进行数值积分。

---

## 5. 物理意义深度解析

### 5.1 矢量效应与标量理论的对比
为了突出矢量的影响，我们对比**标量 Debye 积分公式**：
$$ U(\rho, z) \propto \int_{0}^{\alpha} \underbrace{\sqrt{\cos\theta}}_{P(\theta)} \underbrace{J_0(k\rho\sin\theta)}_{\text{衍射项}} e^{ikz\cos\theta} \sin\theta \, d\theta $$

**对比结论**：
1.  **对称性**：标量公式只包含 $J_0$，意味着焦斑永远是圆对称的（Airy斑）。而矢量公式包含 $I_2 \cos(2\varphi)$ 项，导致 x-线偏振光的焦斑在 x 方向略微拉伸（变成椭圆）。
2.  **纵向分量**：标量理论假设光只有横向振动。矢量理论明确给出了纵向分量 $E_z$（由 $I_1$ 项贡献）。

### 5.2 纵向分量 $E_z$ 与去偏振效应
即使入射光是纯横向偏振，在高 NA 聚焦下也会产生显著的 $E_z$ 分量。
公式项：$E_z \propto -2i I_1 \cos(\varphi)$。

**物理图像分析 (Vector Nodding/Tilting)**：
想象光线从透镜四周汇聚：
*   在 **x-z 平面** ($\varphi = 0^\circ, 180^\circ$)：光线向光轴弯曲（“低头”）。为了保持电场矢量垂直于波矢量，原本沿 x 轴的电场必须向 z 轴倾斜。这产生了一个显著的 $z$ 分量投影。
*   在 **y-z 平面** ($\varphi = 90^\circ, 270^\circ$)：光线虽向光轴弯曲，但电场矢量（沿 x 轴）依然天然垂直于新的传播方向。矢量仅发生了平移，没有发生向 z 轴的倾斜。因此 $z$ 分量为零。

**双瓣结构**：
上述不均匀性（$\cos\varphi$ 依赖性）导致 $E_z$ 呈现**哑铃状（双瓣）结构**。中间被一条暗线（y 轴）分开。

### 5.3 相位异常与古依相移 (Gouy Phase Shift)
积分中的 $e^{ikz\cos\theta}$ 项不仅描述了离焦，还蕴含了古依相移。
当我们观测焦点 ($z=0$) 时，不同 $\theta$ 角的平面波在此处相干叠加。
*   **物理直觉**：当光波被压缩到波长量级的狭小空间（焦斑）时，受到**海森堡不确定性原理**（位置 $\Delta x$ 极小 $\rightarrow$ 动量分布 $\Delta k$ 极大）的约束，光子为了“通过”这个狭隘通道，其相速度（Phase Velocity）会比平面波快。
*   **结果**：这导致在穿过焦点时，波前相对于理想平面波产生了一个额外的 $\pi$ 相位跃变（Phase Anomaly）。这使得几何焦点处的波前并不是平面的，而是存在相位弯曲。

---

## 6. 像差建模与广义光瞳函数

在实际系统中，透镜并不完美，或者样本折射率失配会引入像差。此时圆对称性被打破，$I_0, I_1, I_2$ 简化公式失效。

### 6.1 广义光瞳函数
我们将切趾因子和像差相位合并，定义**广义光瞳函数 $\mathcal{P}(\theta, \phi)$**：
$$ \mathcal{P}(\theta, \phi) = \underbrace{\sqrt{\cos\theta}}_{\text{振幅项}} \cdot \underbrace{\exp[i k \Phi(\theta, \phi)]}_{\text{相位像差项}} $$

其中波像差函数 $\Phi$ 通常使用 **Zernike 多项式** 展开：
$$ \Phi(\rho_{pupil}, \phi) = \sum C_n^m Z_n^m(\rho_{pupil}, \phi) $$
例如：
*   **球差 (Spherical)**: 圆对称，只改变 $I_0, I_1, I_2$ 的数值，不破坏对称性。
*   **像散 (Astigmatism) / 慧差 (Coma)**: 破坏圆对称性，导致 $E_z$ 的双瓣结构扭曲、旋转，不仅影响分辨率，还会严重影响 STED 等依赖零点深度的应用。

### 6.2 数值计算策略：FFT 方法
对于含复杂像差的系统，工程上标准的方法是利用**快速傅里叶变换 (FFT)**。
我们将 Richards-Wolf 积分重写为 2D 逆傅里叶变换形式：

1.  **坐标代换**：将 $(\theta, \phi)$ 映射到空间频率域 $(k_x, k_y)$。
    $k_x = k \sin\theta \cos\phi, \quad k_y = k \sin\theta \sin\phi$
2.  **Jacobian 变换**：$d\Omega \rightarrow \frac{dk_x dk_y}{k k_z}$。
3.  **FFT 公式**：
    $$ \mathbf{E}(x, y, z) = \mathcal{F}^{-1} \left\{ \left[ \frac{\mathbf{a}(k_x, k_y) e^{ik \Phi}}{\cos\theta} \right] e^{ik_z z} \right\} $$

通过在频域网格上计算矢量场 $\mathbf{a}$ 和像差相位，然后执行一次 2D-FFT，即可获得焦平面的完整矢量场分布。

---

## 7. 总结与应用领域

Richards-Wolf 矢量衍射理论是现代纳米光学的基石。
*   **STED/4Pi 显微镜**：需要精确计算甜甜圈光束（Donut beam）中心的剩余强度（由矢量效应和像差引起），这直接决定了超分辨的极限分辨率。
*   **光镊 (Optical Tweezers)**：轴向捕获力（Trapping Force）主要依赖于强梯度的 $E_z$ 分量，标量理论无法预测。
*   **激光微纳加工**：预测透明介质内部加工时的焦斑畸变，指导自适应光学（Adaptive Optics）校正。

该理论完美地将几何光学的偏振追踪与波动光学的干涉积分结合，展示了高 NA 下光的丰富物理行为。

