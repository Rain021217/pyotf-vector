# pyotf-vector 项目详解与完整工作流指南（中文）

> 本文面向三类读者：
> 1) 第一次接触本项目的研究生/工程师；
> 2) 想快速验证模型物理合理性的算法开发者；
> 3) 想扩展共聚焦、矢量激发、4Pi 成像模型的贡献者。

---

## 1. 项目定位与核心能力

`pyotf-vector` 是一个用于**显微成像点扩散函数（PSF）/光学传递函数（OTF）建模**与**相位恢复**的 Python 项目。其核心价值是把经典的 pupil-based FFT 模型与可解释的像差建模（Zernike）结合，并扩展到：

- 宽场（Widefield）
- 共聚焦（Confocal）
- 4Pi-Confocal
- 矢量衍射参考求解（Richards-Wolf / Debye 直接积分，慢但可作为基准）

这使项目同时覆盖：

- 快速工程模拟（FFT 路线）
- 物理一致性对照（矢量 Debye 参考解）
- 像差参数化与逆问题（phase retrieval + Zernike）

---

## 2. 目录结构与模块职责

### 2.1 核心 Python 模块

- `pyotf/otf.py`
  - 定义 `BasePSF`, `HanserPSF`, `SheppardPSF`
  - 支持 `vec_corr`（矢量校正近似）与 imaging condition（sine/herschel）
  - 支持 `apply_aberration` / `apply_named_aberrations`（像差注入）

- `pyotf/microscope.py`
  - 在 PSF 模型层上封装显微镜级模型：
    - `WidefieldMicroscope`
    - `ConfocalMicroscope`
    - `FourPiConfocalMicroscope`
    - `ApotomeMicroscope` / `SIM` 族
  - 共聚焦中提供针孔建模、激发/检测路径分离、像差分路径配置

- `pyotf/vector.py`
  - 提供 **Richards-Wolf / Debye** 参考积分：
    - `richards_wolf_linearly_polarized_psf`
    - `zernike_phase_map`
    - `axial_profile`
    - `optical_sectioning_ratio`
  - 作为 “慢但更接近第一性原理” 的验证通道

- `pyotf/phaseretrieval.py`
  - 迭代相位恢复流程，将实验 PSF 反推 pupil 相位/幅度

- `pyotf/zernike.py`
  - Zernike 多项式与 Noll 序号转换

- `pyotf/utils.py`
  - FFT 包装、坐标工具、数据预处理等通用函数

### 2.2 测试模块

- `tests/otf_test.py`：OTF/PSF 基础性质测试
- `tests/phaseretrieval_test.py`：相位恢复流程测试
- `tests/microscope_test.py`：confocal / 4Pi 行为测试
- `tests/vector_test.py`：矢量参考积分模块测试

### 2.3 教程与文档

- `README.md`：项目简介、安装、模块说明
- `Richards-Wolf_Vector_Diffraction_Summary.md`：矢量衍射理论说明
- `docs/confocal_vector_optimization_plan.md`：confocal/4Pi 优化路线
- `notebooks/`：示例 notebook 与开发记录

---

## 3. 理论到代码的映射

## 3.1 FFT pupil 路线（快）

思想：
1) 构造 pupil（频域）
2) 施加 defocus / condition / vec_corr / aberration
3) IFFT 得到振幅场 `PSFa`
4) 强度 `PSFi = |PSFa|^2`（按分量求和）

优点：快、易扫参、适合工程系统设计。

---

## 3.2 Richards-Wolf / Debye 路线（慢）

`pyotf/vector.py` 通过显式积分实现高 NA 矢量聚焦场计算：

- 在 `(theta, phi)` 角谱上积分
- 采用线偏振入射下矢量分量表达
- 可叠加 Zernike 相位项

优点：物理可解释强，可作为 FFT pupil 近似/实现的对照基准。

---

## 3.3 共聚焦成像模型（当前实现）

在 `ConfocalMicroscope` 中：

- 检测链路：`self.model`（发射波长）
- 激发链路：`self.model_exc`（激发波长）
- 总体 confocal PSF：
  \[
  PSF_{confocal} = PSF_{det,pinhole} \times PSF_{exc}
  \]

其中针孔通过圆孔核建模，可选两种模式：

- `pinhole_mode="object"`：在空间域与检测 PSF 做卷积
- `pinhole_mode="detector"`：在频域对检测链路进行乘法滤波

两者都旨在体现针孔对轴外背景与旁瓣的抑制。

---

## 3.4 4Pi-Confocal（当前实现）

在 `FourPiConfocalMicroscope` 中，激发与检测分别可由双臂场叠加：

- 每条通道支持：
  - 相位控制 `phase_exc` / `phase_det`
  - 振幅比 `amp_ratio_exc` / `amp_ratio_det`
  - 是否相干叠加（interfere flags）
  - arm1/arm2 独立像差

并提供 `phase_scan` 生成相位扫描栈，方便做 I5M/I5S 类后续流程。

---

## 3.5 像差建模

项目支持两条像差输入方式：

1. 系数方式（`mcoefs`, `pcoefs`）
2. 命名方式（`aberrations={"defocus": 0.2, ...}`）

对 confocal / 4Pi 可做到按路径与按臂独立注入，这对于模拟真实系统中激发检测不匹配、双臂不对称非常关键。

---

## 4. 关键类快速索引（建议先读）

### 4.1 `ConfocalMicroscope`

常用参数：

- `wl_exc`：激发波长
- `pinhole_size`：针孔大小（Airy 单位）
- `pinhole_mode`：`object` / `detector`
- `aberrations_em` / `aberrations_exc`
- `pcoefs_em` / `mcoefs_em` / `pcoefs_exc` / `mcoefs_exc`
- `excitation_vec_corr`

常用属性：

- `PSF`：最终归一化系统 PSF
- `excitation_psf`
- `detection_psf`

### 4.2 `FourPiConfocalMicroscope`

新增重点参数：

- `phase_exc`, `phase_det`
- `amp_ratio_exc`, `amp_ratio_det`
- `aberrations_*_arm1`, `aberrations_*_arm2`
- `interfere_excitation`, `interfere_detection`

常用方法：

- `phase_scan(phases, channel="both")`

### 4.3 `richards_wolf_linearly_polarized_psf`

适用于：

- 校验 FFT pupil 路线输出趋势
- 对高 NA 偏振效应进行参考计算
- 验证像差项是否产生预期形变

---

## 5. 建议的新手工作流（从 0 到可复现实验）

## 第一步：安装与测试

```bash
pip install -r requirements.txt
pytest -q
```

如果全部通过，说明本地环境可运行基础建模与测试。

## 第二步：先跑宽场/OTF 直观图

- 使用 `HanserPSF` 生成 PSF、OTF
- 看中心峰、旁瓣、轴向拖尾

## 第三步：跑 confocal 针孔扫描

- 固定光学参数
- 扫 `pinhole_size = [0.2, 0.5, 1.0, 2.0]`
- 观察光学切片能力（离焦/在焦比）

## 第四步：加入像差

- 先加 defocus, astigmatism
- 分别对 excitation 与 detection 加
- 对比不同路径像差导致的 PSF 变化

## 第五步：切到 4Pi

- 扫相位 `phase_exc`/`phase_det`
- 加入 arm1/arm2 像差不平衡
- 观察轴向干涉条纹与峰值移动

## 第六步：用矢量 Debye 参考解做对照

- 选一个小网格，跑 `richards_wolf_linearly_polarized_psf`
- 与 FFT 路线的 profile 做趋势比较

---

## 6. 如何判断“模型是否工作正常”

建议关注如下检查量：

1. `PSF.sum() == 1`（归一化）
2. 减小针孔后离焦背景下降
3. 增大 NA 或启用矢量效应时，横向/轴向结构变化合理
4. 4Pi 相位扫描时，轴向 profile 随相位变化
5. 引入像差后，PSF 不再与无像差相同

---

## 7. 常见问题（FAQ）

### Q1: `pinhole_mode="object"` 和 `"detector"` 为什么结果可能不完全一样？

因为二者是不同实现路径（空间域 vs 频域）的离散近似，受采样与边界条件影响。趋势通常一致，数值可有细微差别。

### Q2: 为什么 `vector.py` 很慢？

它是直接积分，复杂度远高于 FFT pupil。建议仅用于小网格验证，不要直接做大规模参数扫描。

### Q3: 命名像差字符串怎么写？

请使用 `pyotf.otf.name2noll` 支持的名称（例如 `"defocus"`, `"vertical astigmatism"` 等，注意空格）。

---

## 8. 开发建议（后续可继续做）

1. 在 notebook 中加入“FFT 模型 vs Debye 参考”的自动对照图
2. 增加 4Pi phase-step 到 I5M/I5S 合成函数
3. 增加更细粒度 benchmark（时间、内存、误差）
4. 把关键 workflow 做成可复用脚本与 CLI

---

## 9. 快速命令清单

```bash
# 运行所有测试
pytest -q

# 只跑显微镜模型相关
pytest -q tests/microscope_test.py

# 只跑矢量参考求解相关
pytest -q tests/vector_test.py
```

---

## 10. 结语

这个项目最有价值的地方不是“某一个模型公式”，而是：

- **工程速度（FFT）** + **物理参考（Debye）** + **像差参数化（Zernike）** 三者可以闭环；
- 你可以先快速扫设计空间，再用参考积分去校验关键点，最后再接入实验数据做相位恢复。

如果你是初学者，建议按本仓库新增教程 notebook 顺序逐个跑，边改参数边看图，理解会非常快。
