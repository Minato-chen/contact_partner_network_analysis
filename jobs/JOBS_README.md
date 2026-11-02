# 服务器提交

在项目根目录使用pbs.sh或grid_engine.sh，按集群修改队列、内存和时长。环境变量CONTACT_PYTHON指定Python绝对路径，STAGE选择01_center或02_calculate，SIMULATION可选。

先确认01成功再提交02；并行计算结束后单独执行python 04_summary/run.py。模板没有在真实服务器实测，qsub返回编号不等于计算成功。
