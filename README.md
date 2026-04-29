# plantMASST Transcriptomics

**Softwares** 
-	**blast-plus**: v2.16.0: https://github.com/ncbi/blast_plus_docs
-	**orfipy v0.0.4**: https://github.com/urmi-21/orfipy
- **samtools v1.21**: https://github.com/samtools/samtools
- **seqkit v2.3.0**: https://github.com/shenwei356/seqkit/releases
- **sequenceserver v3.1.0**: https://github.com/wurmlab/sequenceserver
- **spades v3.15.5**: https://github.com/ablab/spades
- **sratoolkit v2.10.9**: https://github.com/ncbi/sra-tools
- **trimgalore v0.6.7**: https://github.com/FelixKrueger/TrimGalore
- **EnTAP v2.0.0**: https://gitlab.com/PlantGenomicsLab/EnTAP 
- **Transeq**: https://www.ebi.ac.uk/jdispatcher/st/emboss_transeq (https://doi.org/10.1093/nar/gkae241)

The following scripts were submitted to a high performance computing cluster via SLURM. All softwares above were pre-installed under a Bioinformatics module on the computing cluster except orfipy, which were installed as described below. Sequenceserver was used as a cloud-based version (https://sequenceserver.com/).

# Transcriptome assembly – multiple datasets
**1. Batch SRA-download**
- Configure your cluster for sratools (v2.10.9) commands by running the command before the first software application:
```
./vdb-config -i
```
- Target RNA-seq datasets were primarily paired-end data. For large-scale transcriptome mining, SRA datasets were downloaded in batches of 100 datasets per one directory with the following script. SRA accession numbers of target data were listed in an SRA.txt file.
```
#!/bin/bash
#SBATCH --job-name=sra-download
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./sratools-%j
source /etc/profile.d/http_proxy.sh
module load Bioinformatics
module load sratoolkit/2.10.9-udmejx7
sed -i 's/\r$//' SRA.txt &
cat SRA.txt |parallel -j 10 xargs -n 25 -P 0 fasterq-dump --split-files --outdir /path/to/directory/ & jobs -l
wait
printf "\n...done\n\n"
```

