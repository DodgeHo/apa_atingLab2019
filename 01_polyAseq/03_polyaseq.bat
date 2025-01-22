@echo off
echo "setting python configurations ..."

:: 替换为 E:\apa_atingLab2019\01_polyAseq
set PPD=E:\apa_atingLab2019\01_polyAseq
set PYTHONPATH=%PPD%;%PYTHONPATH%
set PADD_GIT=%PPD%\paddle-git

:: 替换为 F:\ProgramData\Anaconda3\envs\py27
set CONDA_ENV_PATH=F:\ProgramData\Anaconda3\envs\py27
set PATH=%CONDA_ENV_PATH%\bin;%PATH%

echo "running polyAseq pipeline ..."

echo "
python %PPD%\polyAseq_lite.py ^
	-C %PPD%\configs\program.conf ^
	-i %PPD%\01_wkd\fastq ^
	-e DKO ^
	-c HCT ^
	-d %PPD%\01_wkd\comp_group.csv ^
	-o %PPD%\01_wkd\out ^
	-b %PPD%\01_wkd\barcode_list.csv ^
	-p 6 ^
	-r 0
"

python %PPD%\polyAseq_lite.py ^
	-C %PPD%\configs\program.conf ^
	-i %PPD%\01_wkd\fastq ^
	-e DKO ^
	-c HCT ^
	-d %PPD%\01_wkd\comp_group.csv ^
	-o %PPD%\01_wkd\out ^
	-b %PPD%\01_wkd\barcode_list.csv ^
	-p 6 ^
	-r 0
