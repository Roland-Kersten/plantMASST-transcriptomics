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
**2. Batch trimming**
- Trimming can improve precursor peptide assembly and core peptide detection in burpitide biosynthetic studies. RNA-seq data was trimmed with TrimGalore (v0.6.7) default settings in batches of 100 datasets.
- Generate directories for fwd reads and rev reads:
```
mkdir input_data_1
mkdir input_data_2
```
- Move fwd reads to input_data_1/ directory:
```
mv *_1.fastq /path/to/input_data_1/
```
- Move rev reads to input_data_2/ directory:
```
mv *_2.fastq /path/to/input_data_2/
```
- Run the batch TrimGalore-trimming script:
```
#!/bin/bash
#SBATCH --job-name=trimgalore
#SBATCH --account=your_account
#SBATCH --partition=standard 
#SBATCH --array=1-(insert-number-of-datasets-in-input_data_1_directory)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=7g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./trimgalore-%j
module load Bioinformatics
module load trimgalore/0.6.7-ztb2tpz
file1=$(ls ./input_data_1/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
file2=$(ls ./input_data_2/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
trim_galore --cores 4 --paired ./input_data_1/${file1} ./input_data_2/${file2}
```
**3. Batch transcriptome assembly (SPAdes)**
- Generate directories for trimmed fwd reads and trimmed rev reads:
```
mkdir input_data_trimmed_1
mkdir input_data_trimmed_2
```
- Move trimmed fwd reads to input_data_trimmed_1/ directory:
```
mv *_1.fq /path/to/input_data_trimmed_1/
```
- Move trimmed rev reads to input_data_trimmed_2/ directory:
```
mv *_2.fq /path/to/input_data_trimmed_2/
```
- Run SPAdes batch assembly:
```
#!/bin/bash
#SBATCH --job-name=SPAdes
#SBATCH --account=your_account 
#SBATCH --partition=standard
#SBATCH --array=1-(insert-number-of-datasets-in-input_data_trimmed_1_directory)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@mail.com
#SBATCH --mail-type=END
#SBATCH --output=./SPAdes-%j
module load Bioinformatics
module load spades/3.15.5-jhe6qq2
file1=$(ls ./input_data_trimmed_1/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
file2=$(ls ./input_data_trimmed_2/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
spades.py --rna -1 ./input_data_trimmed_1/${file1} -2 ./input_data_trimmed_2/${file2} -o spades_$file1
cd spades_$file1
mv transcripts.fasta /path-to-directory/spades_$file1\.fasta
for f in *_1_val_1.fasta; do
  mv -- "$f" "${f%_1_val_1.fasta}.fasta"
```
**4. BLAST database formatting** 
- BLAST database generation of transcriptome assemblies by Sequenceserver requires reduction of fasta headers to less than 51 letters. The following script was used for SPAdes-specific assembly fasta file formatting for Sequenceserver-based BLAST database generation and commandline PHI-BLAST search.
```
#!/bin/bash
#SBATCH --job-name=rename
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=10g
#SBATCH --mail-user=your@mail.com
#SBATCH --mail-type=END
#SBATCH --output=./rename-%j

for i in *fasta; do n="${i%.fasta}"; sed -i.bak "s/>[^_]\+/>$n/" $i; done
cat *fasta
sed -i 's/_length//g' *fasta
sed -i 's/_cov//g' *fasta
sed -i 's/_g//g' *fasta
sed -i 's/_i//g' *fasta
```
**5. Commandline BLAST search**  
- Commandline PHI-BLAST was applied to searching de novo assembled plant transcriptomes for moroidin core peptides. SPAdes assemblies were combined to a single fasta-file for BLAST database generation prior to PHI-BLAST search.
- a.	Generate transcriptome input file from multiple transcriptome assembly fasta files in a directory:
```
cat *.fasta > all.fasta
```
- b.	Generate query.faa file:
```
nano query.faa
```
- Copy target protein sequence (moroidin cyclase KjaBURP with core peptide sequence QLLVWRGH), for searching moroidin precursor peptides:
```
>QIG55799.1 BURP domain protein [Kerria japonica]
MACRLSLIFAFLCLTLVACHAALSPQEVYWNSVFPQTPMPKTLSALVQPAAKNFIRYKKVDDGQTQDIDV
AADNQLLVWRGHVAIDDDAAADNQLLVWRGHVAIDDDDAAADNQLLVWRGHVAIHDDAAADNQLLVWRAH 
VANDDVDARNLLRKDPSRVLVFLEKDVHPGKTMKYSLIRSLPKSATFLPRNTAESIPFSSSKLLEILIQF 
SVQPKSVEANAMTEAILKCEVPAMRGEAKYCATSLESMIDFVTSRLGRNIRAISTEVEEGATHVQNYTIY 
HGVKKLTDKKVITCHRLRYPYVVFYCHELENTSIYMVPLKGADGTNAKAITTCHEDTSEWDPKSFVLQLL 
KVKPGTDPVCHFLSESDVVWVSNHGTYKPA
```
- c. Generate phi_pattern.txt file:
```
nano phi_pattern.txt
```
- Add the following text for searching a core peptide ‘QLLVWRGH’ and save:
```
PA QQL-x(2)-W
```
- d. PHI-BLAST search
  - Install orfipy before PHI-BLAST search as described above. The same orfipy translation parameters were applied for PHI-BLAST search as for Sequenceserver-based tblastn search of burpitide cyclase sequences (i.e. minimum of 450 bp open reading frame length, translation between stop codons).
```
#!/bin/bash
#SBATCH --job-name=phiblast-QLLVWRGH
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@mail.com
#SBATCH --mail-type=END
#SBATCH --output=./phiblast-%j
module load Bioinformatics
module load blast-plus/2.16.0
module load seqkit2
orfipy --pep all.pep --min 300 --between-stops all.fasta
cp orfipy_all.fasta_out/all.pep .
awk '/^>/ {print $1; next} {print}' all.pep > cleaned_all.pep
makeblastdb -in cleaned_all.pep -dbtype prot -out cleaned_all_db
psiblast -db cleaned_all_db -query query.faa -out cleaned_all_phiblast.txt -phi_pattern phi_pattern.txt
grep ">" cleaned_all_phiblast.txt | awk '{print $1}' | tr -d '>' > hits.txt
seqkit grep -f hits.txt cleaned_all.pep > phiblast_all_sequences.pep
```
