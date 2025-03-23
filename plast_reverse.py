import sys
import subprocess
import pkg_resources

# 检查并安装 Biopython
required_packages = {'biopython'}
installed_packages = {pkg.key for pkg in pkg_resources.working_set}
missing_packages = required_packages - installed_packages

if missing_packages:
    print(f"Missing packages: {missing_packages}, installing now...")
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', *missing_packages])

# 导入 Biopython
from Bio import AlignIO

# 其余代码
print("Biopython successfully imported!")



# 读取比对结果
from Bio import SeqIO
from Bio.Seq import reverse_complement

# 读取叶绿体基因组
genome = SeqIO.read("embplant_pt.K105.complete.graph1.1.path_sequencewhite.fasta", "fasta")

# SSC 区域起止位置（Python 索引为 0 基准）
SSC_start, SSC_end = 114240, 132862  # 你的 SSC 区域

# 提取 SSC 区域序列并反转互补
SSC_original = genome.seq[SSC_start:SSC_end]
SSC_reversed = reverse_complement(SSC_original)

# 替换原 SSC 区域
genome.seq = genome.seq[:SSC_start] + SSC_reversed + genome.seq[SSC_end:]

# 保存修改后的基因组
SeqIO.write(genome, "chloroplast_genome_SSC_corrected.fasta", "fasta")

print(f"SSC 区域 ({SSC_start+1}-{SSC_end}) 已成功反转！")

