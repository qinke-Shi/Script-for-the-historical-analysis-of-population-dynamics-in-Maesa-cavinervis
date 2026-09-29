# 切换到工作目录
cd /data_forrest/xubo/shiqinke/Maesa_pop/compare_gene/2_mcmTree

# 设置软件路径
MCMCTREE=/home/xubo/Software/paml4.9j/src/mcmctree

# 打印开始时间
echo "Job started at: $(date)"

# 运行 MCMctree
$MCMCTREE mcmctree.ctl

# 打印结束时间
echo "Job finished at: $(date)"
