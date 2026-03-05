# Confocal / Vector-PSF / 4Pi 优化计划（结合 Richards-Wolf 理论）

## 1. 当前模型物理一致性检查结论

### 1.1 Confocal 三维成像理论一致性
- 现有 `ConfocalMicroscope` 使用 `PSF_confocal = PSF_det,pinhole ⨉ PSF_exc` 的乘积模型，符合经典共聚焦在荧光成像（照明与检测独立）下的一阶近似。
- 针孔通过对检测 PSF 在横向与圆孔核卷积实现，物理上等价于对探测面有限孔径的积分；孔径减小时，轴外背景与旁瓣抑制增强。
- 局限：原实现仅做“检测路径卷积”，未显式支持矢量激发，也未提供检测/激发独立像差接口。

### 1.2 Zernike 三维像差建模合理性
- 项目已有 `apply_aberration`（Hanser pupil + Zernike）属于广义光瞳函数方法，形式与 Richards-Wolf 文档中
  `P(θ,φ)=A(θ,φ)exp(i k Φ(θ,φ))` 一致。
- 优点：可控、可解释、可与相位恢复结果闭环。
- 局限：当前像差入口主要聚焦在单模型对象，缺少显微成像级（exc / det / 4Pi 双臂）统一接口。

### 1.3 dual 机制与 4Pi 可扩展性
- `SheppardPSF(dual=True)` 通过双壳构造体现“对向物镜采集几何”。
- 但若要物理复现 4Pi，必须显式处理**双臂复振幅相干叠加**与相位控制（`Δφ`）。
- 因此建议在显微镜层新增 4Pi-confocal 模型：分别对 excitation / detection 计算前后向场并干涉，再与针孔联动。

## 2. 已完成的第一阶段改造

1. 为 `ConfocalMicroscope` 增加了激发/检测像差参数接口：
   - `pcoefs_em`, `mcoefs_em`, `pcoefs_exc`, `mcoefs_exc`
2. 增加矢量激发控制参数：
   - `excitation_vec_corr`（默认 `"total"`），并提供 `excitation_psf` 接口。
3. 新增 `FourPiConfocalMicroscope`：
   - 支持 excitation / detection 双通道相位参数 `phase_exc`, `phase_det`
   - 支持是否干涉开关 `interfere_excitation`, `interfere_detection`
   - 在检测路径继续应用针孔卷积，实现 4Pi-confocal 的基础模型闭环。

## 3. 分阶段优化路线图（建议按顺序执行）

### 阶段 A：理论与数值一致性验证
1. 建立 benchmark：
   - widefield / confocal / 4Pi-confocal 的轴向与横向剖面对比
2. 验证 pinhole 作用：
   - 扫描 0.2–2 AU，量化旁瓣积分、离焦平面能量、FWHM
3. 验证矢量激发：
   - `vec_corr=none` vs `total`，比较高 NA 下焦斑偏振非对称与 `Ez` 贡献

### 阶段 B：像差框架统一
1. 统一 API：
   - `aberration_exc`, `aberration_det`, `aberration_arm1`, `aberration_arm2`
2. 支持双臂独立 Zernike：
   - 模拟 4Pi 臂间像差不匹配导致的条纹对比下降与焦点漂移
3. 增加像差组合工具：
   - named-aberration 到系统级参数映射

### 阶段 C：4Pi 物理完备化
1. 双臂幅度不平衡（振幅比）
2. 臂间群时延/色散近似项
3. 相位扫描模式：生成 I5M / I5S 风格数据接口

### 阶段 D：性能优化
1. 复用 pupil / defocus 缓存（避免重复 FFT）
2. 对 pinhole 卷积采用 OTF 域乘法可选路径
3. 对参数扫描任务增加批量接口（vectorized）

### 阶段 E：测试与文档
1. 增加单元测试：
   - 针孔收缩导致离焦背景下降
   - 4Pi 相干开关与相位调节生效
2. 增加示例 notebook：
   - `Vector Excitation + 3D Confocal`
   - `4Pi Confocal + Aberrations`
3. 输出“理论-代码映射表”：逐条对应 Richards-Wolf 核心公式与代码位置

## 4. 后续可选增强
- 引入严格矢量 Debye 积分（θ, φ 直接积分）作为高精度参考解；当前 FFT pupil 法保留为快速工程实现。
- 在 confocal 检测端加入“针孔共轭面采样”模型（而非仅物空间卷积近似）用于大针孔时的误差校正。
