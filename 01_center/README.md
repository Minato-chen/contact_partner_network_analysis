# 轨迹居中

只负责从原始模拟中选帧、把slab移到盒中心并恢复完整链，不承担原始轨迹平衡或质量审查。原始文件保持只读。

- run.py：批量入口与参数；dependencies/center_local.py：居中算法。
- 输入：simulations.yaml中的原始拓扑和轨迹。每segment一条有序等长链，每残基一珠，固定正交盒。
- 方法：先按密度峰平移，再按最大连续占据段的加权中心调整，最后恢复整条链。
- 默认参数：原始零起始帧号start=750、step=3、stop=None，预期750帧；在run.py顶部修改。

```bash
python -u 01_center/run.py
python -u 01_center/run.py --simulation B01_T260
```

结果位于output/<模拟编号>/：centered.dcd为居中轨迹，topology.pdb和chains.json描述拓扑与链，frame_index.csv对应原始帧号，centered.diagnostics.csv保存居中诊断，metadata.json记录来源和参数。全部完成后才有COMPLETE。

01保存坐标单位为Å，后续MDTraj读取时转换为nm。若清单提供原始保存帧间隔time_step_ps，可计算实际时间；缺少则时间留空，不猜测。上游更改后需要重算后续结果。
