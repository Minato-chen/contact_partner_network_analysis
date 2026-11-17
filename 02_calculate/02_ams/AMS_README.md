# $A_{\mathrm{MS}}$ 序列谱分解与分子内外接触耦合分析 (`02_ams`)

本模块作为 `02_calculate` 下的第二个核心计算任务，专门用于连接**一级序列电荷矩阵特征 ($A_{\mathrm{MS}}$)** 与**三维分子动力学轨迹中的分子内/分子间接触图谱 (`01_contacts`)**，为论文第四部分提供完整的定量因果证据链。

---

## 一、 模块架构与输出层级

```text
02_calculate/02_ams/
├── run.py                 # 自动调度 01_spectrum 与 02_coupling
├── 01_spectrum/           # 序列级：SCD 符号拆解、多尺度谱特征与残基主导模式权重
│   ├── run.py
│   └── spectrum.py
├── 02_coupling/           # 耦合级：将序列电荷主导模式与 01_intra / 02_inter / 03_joint 对齐
│   ├── run.py
│   └── coupling.py
└── AMS_README.md
```

结果输出至 `03_results/<simulation_id>/03_ams/`，经 `04_summary/run.py` 自动聚合成三张横向对比宽表：
- `04_summary/output/ams.csv`（来自 `01_spectrum` 的纯序列分解指标）
- `04_summary/output/ams_dense_chains.csv`（浓相系综序列-接触耦合指标）
- `04_summary/output/ams_central_chain.csv`（中心探针链序列-接触耦合指标）

---

## 二、 子模块 1：序列电荷谱分解 (`01_spectrum`)

1. **SCD 符号拆解 (`scd_decomposition.csv`)**：
   - `spectrum.scd_pos_pos` ($\mathrm{SCD}_{++} > 0$)：正-正同号电荷对贡献；
   - `spectrum.scd_neg_neg` ($\mathrm{SCD}_{--} > 0$)：负-负同号电荷对贡献；
   - `spectrum.scd_same_sign` ($\mathrm{SCD}_{++} + \mathrm{SCD}_{--}$)：同号排斥总贡献；
   - `spectrum.scd_pos_neg` ($\mathrm{SCD}_{+-} < 0$)：正-负异号吸引贡献。
2. **多尺度电荷矩阵谱特征 (`spectrum_by_scale.csv`, $\ell \in \{2, 5, 10, 20\}$)**：
   - `spectrum.f_pos_ms`, `spectrum.f_neg_ms`：$\ell \in \{5, 10, 20\}$ 对数积分的正、负电荷谱集中度；
   - `spectrum.ams` ($A_{\mathrm{MS}} = (\Delta(5) + 2\Delta(10) + \Delta(20))/4$)：多尺度正负谱集中度差；
   - `spectrum.cms` ($C_{\mathrm{MS}} = (C(5) + 2C(10) + C(20))/4$)：多尺度总体谱集中度；
   - `spectrum.pos_principal_center`, `spectrum.neg_principal_center`, `spectrum.cross_principal_center_separation`（$\ell=10$ 主导块序列位置与分离度）。
3. **残基主导模式权重 (`residue_modes.csv`)**：
   - 每个残基位点 $i \in [0, 128]$ 在各尺度下的主特征向量权重 $w_+(i) = v_{+,1}^2(i)$ 与 $w_-(i) = v_{-,1}^2(i)$。

---

## 三、 子模块 2：序列-接触联合耦合分析 (`02_coupling`)

将 `01_spectrum` 与 `01_contacts`（`01_intra`、`02_inter`、`03_joint`）严格配对，输出：

1. **分子内接触的电荷符号与距离拆解 (`intra_charge_decomposition.csv`)**：
   - `coupling.intra_attract_contacts` (`pairs/chain`)：链内异号电荷（`+-`）自中和接触对数；
   - `coupling.intra_long_range_attract_contacts` (`pairs/chain`)：链内长程（$|i-j| \ge 12$）异号自锁接触对数；
   - `coupling.intra_attract_fraction` (`dimensionless`)：链内异号接触占总链内接触的比例；
   - `coupling.inter_attract_partition` (`dimensionless`)：$\frac{\text{inter.attract\_contacts}}{2 \times \text{intra\_attract\_contacts} + \text{inter.attract\_contacts}}$，量化异号静电吸引投向“分子间交联”而非“分子内自中和”的分配比例；
   - `coupling.pos_mode_intra_degree`, `coupling.neg_mode_intra_degree`：正、负主导电荷块的加权链内接触度。
2. **分子间多价转化与单伙伴垄断度 (`summary.csv` & `frame_means.csv`，带 50 帧块 SEM)**：
   - `coupling.max_partner_attract_share [threshold=k]`：$\frac{\text{max\_attract\_contacts\_per\_partner}}{\text{attract\_contacts}}$，单一最强伙伴垄断全链异号接触的份额（严格 $\in [0, 1]$，直接表征静电拉链锁死程度）；
   - `coupling.partner_per_attract_contact [threshold=k]`：$\frac{\text{partner\_count}}{\text{attract\_contacts}}$，每对分子间异号接触换来的有效伙伴数（多价转化效率）；
   - `coupling.partner_per_neg_residue [threshold=k]`：$\frac{\text{partner\_count}}{N_-}$，每个少数派负电荷（限制性贴纸）平均招募的有效伙伴链数。
