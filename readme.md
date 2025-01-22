# Download APA reproducible package

```
mkdir -p $HOME/projects
cd $HOME/projects
git clone https://github.com/DodgeHo/apa_atingLab2019.git
```

- In this guideline, the package base directory is "$HOME/projects/apa_atingLab2019". In the following steps, it is necessary to modify scripts to be matched with your local setup.
- 注意apa_atingLab2019最好安装在 Linux下的~/projects/apa_atingLab2019路径，如果不是，下文有些地方请注意修改

# Prerequisite

1. python 3 与 conda环境：
   
   由于python2.7很多包已经不再适用，建议使用python3
   
   ```
   # 如果电脑已有conda，跳过；如果没有conda，先安装一个miniconda;
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   sudo chmod +x Miniconda3-latest-Linux-x86_64.sh
   sudo ./Miniconda3-latest-Linux-x86_64.sh -b -f -p /usr/local
   #请注意sudo可能需要输入密码
   
   # 如果电脑已有pip，跳过；如果没有pip，需要安装
   sudo apt install python3
   sudo apt install python3-pip
   ```

2. conda环境设置，并安装conda和pip的扩展包：
   
   ```
   # 进入代码文件夹
   cd apa_atingLab2019
   
   # 创建 Conda 环境（例如命名为 bioinfo，以后只要activate就可以用了）
   # 这一步是为了避免一些已有的安装冲突，因此很重要
   conda create -n bioinfo  -y
   conda activate bioinfo
   
   # 安装pip的扩展包
   pip install python_requirements.txt 
   # 安装conda的扩展包
   programs=(
       "mosdepth" "sra-tools" "star" "xz" "wget"
       "configparser" "numpy" "pandas" "pysam" "argparse"
       "pyfaidx" "Bio" "deeptools" "regex"
   )
   for prog in "${programs[@]}"; do
       conda install -c bioconda "$prog" -y
   done
   ```

3. 安装R语言和R的扩展包
   
   ```
   conda install -c bioconda r-base=3.4.3 -y
   
   # 定义要安装的 R 包列表
   programs=(
       "goldmine" "argparse" "BSgenome.Hsapiens.UCSC.hg19" "corrplot"
       "data.table" "DEXSeq" "e1071" "edgeR" "ggbio" "ggplot2" "gridExtra"
       "IlluminaHumanMethylation450kanno.ilmn12.hg19" "matrixStats" "openxlsx"
       "parallel" "reshape" "RPMM" "Rsamtools" "Rsubread" "stringr" "ggseqlogo"
       "JASPAR2018" "TFBSTools" "UpSetR" "VennDiagram" "wateRmelon" "devtools"
   )
   
   # 遍历程序列表并安装
   for prog in "${programs[@]}"; do
       Rscript -e "install.packages('$prog', repos='https://cloud.r-project.org/')"
   done
   
   # 克隆 Handy仓库并安装 （handy这个包是非官方的，所以需要自己下载安装）
   cd ..
   !git clone https://github.com/jeffbhasin/handy.git
   cd handy
   !Rscript -e "library(devtools); install_github('jeffbhasin/handy')"
   cd ..
   ```

4. 安装其他需要的程序
   
   ```
    安装 bowtie2,samtools,Java（picard 依赖）
   apt-get update
   apt-get install -y bowtie2
   apt-get install -y samtools
   apt-get install -y openjdk-11-jre
   
   # 下载 picard
   wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar
   
   # 下载 UCSC 工具集
   wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
   wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
   
   sudo chmod +x twoBitToFa
   sudo chmod +x wigToBigWig
   ```

5. 将R语言需要的项目工作环境加到环境变量里。请注意，如果你的项目路径不是~/projects/apa_atingLab2019路径，请在下面修改你的路径到环境变量里。
   
   ```
   # 设置环境变量
   export R_UTIL_APA="$HOME/projects/apa_atingLab2019/resource/libs"
   export PADD_GIT="$HOME/projects/apa_atingLab2019/01_polyAseq/paddle-git"
   # 将环境变量写入 ~/.bashrc（可选）
   echo "export R_UTIL_APA=$R_UTIL_APA" >> ~/.bashrc
   echo "export PADD_GIT=$PADD_GIT" >> ~/.bashrc
   # 使 ~/.bashrc 生效
   source ~/.bashrc
   ```

# PolyA-seq data processing

接下来我们做运行的准备工作

6. 打开配置文件 (`apa_atingLab2019/01_polyAseq/configs/program.conf`) 并编辑项目路径，需要修改的路径有五个，把~/projects/apa_atingLab2019修改成你的项目路径。（如果你是按照说明要求安装的位置，则不用修改）

```bash
[bowtie2-build]
bin_path = bowtie2-build
ref_hg19_path = ~/projects/apa_atingLab2019/resource/ref/hg19_PhiX.fa
ref_hg19_prefix = ~/projects/apa_atingLab2019/resource/ref/hg19_PhiX
chrom_size = ~/projects/apa_atingLab2019/resource/ref/ChromSizes.hg19.txt

[ensembl]
biotype_table = ~/projects/apa_atingLab2019/resource/ensembl/ensembl_gene_biotypes.csv
gtf_file = ~/projects/apa_atingLab2019/resource/ensembl/Homo_sapiens.GRCh37.87.gtf
```

7. 修改脚本 (`03_polyaseq.sh`,位于`apa_atingLab2019/01_polyAseq/03_polyaseq.sh`) 的路径，规则同上，需要修改的只有一行
   
   ```
   export PPD=$HOME/projects/apa_atingLab2019/01_polyAseq
   ```

8. 现在我们可以开始运行脚本了
   
   ```bash
   # 先确保自己位于子文件夹01_polyAseq
   cd ~/projects/apa_atingLab2019/01_polyAseq
   bash ./01_resource_prep.sh
   ```

9. 运行脚本2下载 PolyA-seq FASTQ 文件.  
   
   ```
   bash ./02_get_fastq.sh
   ```
   
   Note that, as of 05/21/2019, the SRAs, 
- SRP083252 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86178)

- SRP083254 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86180) 
  
  is not public. Secure tokens are available only reviewers to access processed data. Any raw data is not available yet until it becomes public.
  
  如果不希望全部下载，可以打开02_get_fastq.sh修改。在其中找到一行
  
  ```
   #srrs=(SRR4091084 SRR4091085 SRR4091086 SRR4091087 SRR4091088 SRR4091089 SRR4091090 SRR4091091 SRR4091104 SRR4091105 SRR4091106 SRR4091107 SRR4091108 SRR4091109 SRR4091110 SRR4091111 SRR4091113 SRR4091115 SRR4091117 SRR4091119)
    #修改为
    srrs=(SRR4091084)
  ```
10. 现在，可以执行polyA-seq processing pipeline（脚本3）了
    
    ```
    bash ./03_polyaseq.sh
    ```

11. Transcription factor binding analysis in the genomic regions in between polyA sites at each 546 APA gene
    
    ```
    Rscript ./04_enrich_perms_100k_select.r
    ```

```
11. Plot sequence logos for the highly enriched TFs
```

 Rscript ./05_plot_seqlogos.r

```
# RNA-seq data processing

1. Contact authors to download RNA-seq FASTQ files.
```

   cd apa_atingLab2019/02_mRNAseq
   bash ./01_get_fastq.sh

```
2. Run STAR 
```

   bash ./02_run_star.sh

```
3. Collect a basic alignment statistics 
```

   Rscript ./03_get_alignment_stats.r

```
4. Analyzing the APA regulatory gene expressions. If the number of CPU's available is 8, 
```

   bash ./04_apafactorexp.sh 8

```
# ChIP-seq data processing

1. We use ENCODE TF and Histone ChIP-Seq processing pipeine from Kundaje lab (e.g., chipseq.bds.20180726_175025_830). If you don't have the pipeline, then, 

1. Download the pipeline at https://github.com/kundajelab/chipseq_pipeline
2. Follow the installation instruction. The installation directories used in this guideline are,
   - `$HOME/apps/chipseq/chip-seq-pipeline2` # ENCODE pipeline path
   - `$HOME/.bds` # bds path
   - `$HOME/apps/miniconda3/bin` #miniconda3 path
   - `/opt/bds_pipeline_genome_data/aquas_chipseq_species.conf` #species configuration file
   - Make sure that the installation is successful
3. Come back to the APA project directory, change to ChIP-seq working directory, and download ChIP/MBD-seq FASTQ files.
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131606

- https://www.ncbi.nlm.nih.gov/sra/SRP001414

  As of 05/21/2019, note that ChIP-seq raw fastq files are not in public.
```

  cd $HOME/projects/apa_atingLab2019/03_chipseq
  bash ./00_get_fastq.sh

```
1. Open the bash script (`chipseq_pipeline_arg.sh`) to modify the paths matched with your local setup.
```

   progd="$HOME/apps/chipseqs/TF_chipseq_pipeline"
   species_conf="/opt/bds_pipeline_genome_data/aquas_chipseq_species.conf"

   export PATH=$HOME/.local/bin
   export PATH=$HOME/.bds:$PATH
   export PATH=$HOME/apps/miniconda3/bin:$PATH
   export PATH=/usr/lib64/ccache:/usr/local/bin:/usr/local/sbin:/usr/bin:/usr/sbin:$PATH

```
2. Launch the ENCODE ChIP-seq pipeline. Note that it will take for a while to complete the job.
```

   bash ./01_chipseq_by_encode.sh 2>&1 | tee logs/01_chipseq_by_encode_$(date '+%Y-%m-%d-%H').log

```
3. Process MBD-seq FASTQ files.
```

   bash ./02_mbdseq_by_bt2.sh 2>&1 | tee logs/02_mbdseq_by_bt2_$(date '+%Y-%m-%d-%H').log

```
4. Generate bigWig files for ChIP-seq alignment files for visualization. If you have a ftp server to host bigWig files to import to UCSC genomebrowser, open the script and modify `ftp://account:passwd@hostname/apa_atingLab2019/chipseq_bigwigs`.
```

   bash ./03_chipseq_tags_to_bigwig.sh
   bash ./04_input_tags_to_bigwig.sh

```
5. Collect some basic alignment statistics
```

   bash ./05_get_align_stats.sh

```
6. Perform a differential ChIP-seq binding site analysis using MANorm.
```

   bash ./06_manorm_diff_chipseq.sh
   Rscript ./06a_merge_manorm_report.r
   Rscript ./06b_store_manorm_sites.r

```
7. Visualization with ChIP-seq differential binding sites. Note that b03_nmf.sh requires that MATLAB is installed and should be in $PATH.
```

   Rscript ./b01_get_apa_intby_chipseq.r
   Rscript ./b02_uniq_xcvg_reg.r
   bash ./b03_nmf.sh
   Rscript ./b04_gen_heatmap_filtw.r
   Rscript ./b05_cluster_analysis.r

```
# TCGA

TCGA data access permission is required to complete this analysis. Consider the following scripts to figure out which procedures/parameters were used in the paper for reference. Contact us for more help.

1. Change to TCGA working directory and download a preprocessed data (TCGA methylation beta values and pA usage ratio on APA genes of interest).
```

   cd $HOME/projects/apa_atingLab2019/04_tcga
   bash ./01_resource_prep.sh
   less -S 546goi_all/github_Supplementary_Table_S6.tsv #table

```
2. Access TCGA RNA-Seq BAM files and 4500 Infinium methylation array. Compute predicted polyA usage and normalized methylation level (beta value). Visualize the correlation with adjusted P-value in scatter plots and gene track along with ChIP-seq binding site read pileup from the cell line model designed.
```

   ./02a_tcga_analy_apa_paur.r
   ./02b_tcga_analy_apa_paur_plot.sh

```
3. To summarize/visualize the correlation between polyA usage predicted from mRNA-seq and methylation level,
```

   cd ./03_tcga_summary
   Rscript ./01_collapse_corr_mpts_max_corr_per_cohort.r
   Rscrtip ./02_summary_corr_table.r

```
4. To generate a box plot of polyA usage and methylation level in HEATR2 per tumor stage.
```

   cd ../04_tcga_clinical_info
   bash ./runme.sh

```
# Reference
```
