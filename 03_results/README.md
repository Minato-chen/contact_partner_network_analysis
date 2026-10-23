# 数值结果

第一层为模拟编号，内部01_basic、02_contacts等分类对应02_calculate。详细计算成功后有COMPLETE和metadata.json；common子目录是对应计算保存的底层数据，不是独立任务。

最终指标另存scalars.json供04_summary收集；单值计算可直接保存<指标>.scalars.json，无需额外目录。详细数据保留在这里，汇总表在04_summary/output，Notebook和图片在05_plot。

模拟结果不进入Git，请另行备份。不要直接修改数值或在这里混入其他模拟的数据。
