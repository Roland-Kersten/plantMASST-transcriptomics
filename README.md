# plantMASST Transcriptomics

**Softwares** 
-	**blast-plus**: v2.16.0: https://github.com/ncbi/blast_plus_docs
-	**orfipy v0.0.4**: https://github.com/urmi-21/orfipy
- **samtools v1.21**: https://github.com/samtools/samtools
- **seqkit v2.3.0**: https://github.com/shenwei356/seqkit/releases
- **sequenceserver v3.1.0**: https://github.com/wurmlab/sequenceserver
- **spades v3.15.5**: https://github.com/ablab/spades
- **sratoolkit v2.10.9**: https://github.com/ncbi/sra-tools
- **kallisto v0.46.0**: https://github.com/pachterlab/kallisto
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
a.	Generate transcriptome input file from multiple transcriptome assembly fasta files in a directory:
```
cat *.fasta > all.fasta
```
b.	Generate query.faa file:
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
c. Generate phi_pattern.txt file:
```
nano phi_pattern.txt
```
- Add the following text for searching a core peptide ‘QLLVWRGH’ and save:
```
PA QQL-x(2)-W
```
d. PHI-BLAST search
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
# Differential gene expression analysis
1.	Batch SRA-download
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
with SRA.txt:
```
nano SRA.txt
```
including sequences:
```
SRR17400452
SRR17400460
SRR17400461
SRR17400462
SRR17400463
SRR17400464
SRR17400465
SRR17400466
SRR17400467
```
2. Batch trimming
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
Run batch TrimGalore-trimming script:
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
3. Combined transcriptome assembly (SPAdes)
- Combine forward reads of all trimmed paired datasets.
```
cat *_1.fq > Aglaonema_1.fq
```
- Combine reverse reads of all trimmed paired datasets.
```
cat *_2.fq > Aglaonema_2.fq
```
- Run SPAdes assembly of trimmed combined fastq files:
```
#!/bin/bash
#SBATCH --job-name=SPAdes
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH --mem=180g
#SBATCH --mail-user=your@mail.com
#SBATCH --mail-type=END
#SBATCH --output=./SPAdes-%j
module load Bioinformatics
module load spades/3.15.5-jhe6qq2
spades.py --rna -1 ./Aglaonema_1.fq -2 ./Aglaonema_2.fq -o spades_Aglaonema
cd spades_Aglaonema
mv Aglaonema.fasta /path/to/directory/Aglaonema.fasta
```
4. BLAST search of candidate burpitide precursors
- BLAST database formatting - BLAST database generation of transcriptome assemblies by Sequenceserver requires reduction of fasta headers to less than 51 letters. Below is a script for SPAdes assembly fasta file formatting for Sequenceserver-based BLAST database generation.
```
#!/bin/bash
#SBATCH --job-name=rename
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
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
- For burpitide precursor search, formatted Aglaonema.fasta file was uploaded to a cloud-based Sequenceserver (v3.1.0) and searched by tblastn (BLAST+ v2.16.0) for homologs of SkrBURP (fused cyclopeptide alkaloid precursor peptide from Selaginella kraussiana) with the following BLAST parameters: evalue 1e-05, matrix BLOSUM62, gap-open 11, gap-extend 1, filter L, max_target_seqs 500.
```
>QXY82431.1 BURP [Selaginella kraussiana]
MAQSLILLFLIVMGCYDVGAIERVHSAKEGVSEAMNDQGHPTIANKAVPASVEDPSQADADNILLYPSYA
AKKAVPEDPSQADADNVLFYRSYAAKKAVPASVEDPSHDADNVLFYPSYVAKKAVPASVEDPSQADADNF
LLYPYAAKKAVPASVEDPSHDADNVLFYPSYVAKKAVPASVEDPSQADADNFLFYPYAAKKAVPASVEDP
SHDADNVLFYPSYVAKKAVPASVEDPSHDADNVLFYPSYVAKKAVPASVEDPSQADADNVLFYPSYAARP
KAPSTMDHDATHNAAQNMHHQKGALLFFRMNTLKVGGEVVIPSLASHALGGRKLMTPRLEAMLGSMELPR
LLSVLKIPLDSTLARKSKTNLGNCRSPPLQGEKKACVSSIRSMTNFARSVLEDKKPLEHLNPSRPVSSPE
HFKIMDVKVVTENSVVTCHPMVFPYALYMCHFVPKSVPIKVTLQDDHDKLVVVPVMCHMDTSEFDPSHLS
FKILNTKPGEAEMCHWMPNSHIMWYTSDGKTRDVL
```
- Download “FASTA of all hits”, 6frame translate the hit transcript sequences with orfipy or with EMBL TRANSEQ .
- AcoBURP-FLLY transcript was identified by searching for core peptide FLLY among translated sequences which led to identification of following transcript-derived BURP-domain-containing protein sequence:
```
>Aglaonema_40161_2574_2845.588964153010_1
VDFATKPFLLYSKGGDGVDHNKMNAGVVSGVDGKPVTVDFSTKPFLLYSKGGNGVDHNKL
NAGVVSGVDGKPVTVDFATKPFLLYSKGGDGVDHNKMNAGVVSSVDGKPVTIDFSTKPIL
LYNKGGDGVDHNKLNAGVVSGVDGKPVTVDFATKPILLYSKGGDGVDHNKLNAGVVSGVD
GKPVTVDFSTKPILLYSKGGDGVDHNKLNAGVISGVNGNPVTVDFATKPFLLYSKGGDGV
DRNKMNAGVVSGVDGKPITVDFSTKPFLLYNKGGDGVDRTKLNAGVVSGVDGKSVTVDFA
TKPFLLYNKGGDGVDHNKMNAGVVSGVDGKHVTVDFSTKPSLLYNKGGDGVDHNKLNAGV
VSGVDGKPVTVDFATKPFLLYSKGGDGVDHNKMNAGVVSGVDGKPVTVDFSTKQILLYNN
GGNGVDHNKLNAGVVSGVDGKPVTVDLATKPFLLYSKGGDGVEHNKLNAGVISRVDGKLV
YSRKASIGSRLSLGGELAAGRNSITHDNIGNPKALNFFLEKDLHVGAKMKMTLDRGMSGD
TFLHREIADTIPFSFNKFPEILRRFAIQPNSPMAEKLEDTLHKCEDLPAASEKKFCATSL
ESMVDFITSSMGTNDLETLETVVNRDAPEQLYTVQSWVVVPTPNWTATICHVYKYPYTVF
HCHWIPQSKPYILKTFREDGSMMDVVAICHVDTSEWPPTFHALQAFGLRPGEATLCHLMS
EGDIVWFPRGARATVTAM
```
5. Kallisto read quantification
-	Index assembled de novo transcriptome with the following kallisto script:
```
#!/bin/bash
#SBATCH --job-name=kallisto_index
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@mail.com
#SBATCH --mail-type=END
#SBATCH --output=./kallisto-index-%j
module load Bioinformatics
module load kallisto/0.46.0
kallisto index -i Aglaonema.idx Aglaonema.fasta
```
- Quantify reads of each paired trimmed fastq dataset (SRA#) against indexed assembled transcriptome with the following kallisto script:
```
#!/bin/bash
#SBATCH --job-name=kallisto_quant_SRA#
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@mail.com
#SBATCH --mail-type=END
#SBATCH --output=./kallisto-quant-%j
module load Bioinformatics
module load kallisto/0.46.0
kallisto quant -i Aglaonema.idx -o output_SRA# SRA#_1_val_1.fq SRA#_2_val_2.fq
```
- Open resulting tsv-files in Microsoft Excel and sort the column ‘target_id’ from ‘AZ’, i.e. transcript names in each SRA#-specific tsv file.
- Generate a new excel sheet with the sorted ‘target_id’ column and all ‘tpm’ columns from each SRA#.tsv file. Add the SRA# to the top of each corresponding tpm column in the newly combined sheet.
- Copy row of AcoBURP-FLLY transcript to the top of the table as the bait gene for correlation analysis.
- Calculate Pearson correlation coefficients of all rows below AcoBURP-FLLY.
- Sort all Pearson correlation coefficient values from highest to lowest.
- Copy all “target_id” transcript IDs with Pearson coefficient values 1-0.7 to a txt file (pearson07_hits.txt) for subsequent annotation on a computational cluster with EnTAP (v2.0.0).
6. Gene annotation with EnTAP
- Compile sequences of top correlation transcripts from txt file and transcriptome fasta file with seqkit2 script:
```
#!/bin/bash
#SBATCH --job-name=seqkit2
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=14g
#SBATCH --mail-user=your@mail.com
#SBATCH --mail-type=END
#SBATCH --output=./seqkit-%j
module load Bioinformatics
module load seqkit/2.3.0-hsk744b
seqkit grep -r -f pearson07_hits.txt Aglaonema.fasta > Aglaonema_pearson07.fasta
```
- Annotate the sequences with positive correlation values (1-0.7) in new fasta file (Aglaonema_pearson07.fasta) with EnTAP. EnTAP requires the installation and configuration of a reference database such as Uniprot. Please refer to EnTAP gitlab (https://gitlab.com/PlantGenomicsLab/EnTAP) for database installation.
```
#!/bin/bash
#SBATCH --job-name=entap
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=14g
#SBATCH --mail-user=your@mail.com
#SBATCH --mail-type=END
#SBATCH --output=./entap-%j
module load Bioinformatics
module load EnTAP/2.0.0
EnTAP --runN --state 4x --run-ini /scratch/rkersten_root/rkersten0/rkersten/Aglaonema_red/entap_run.params --entap-ini /scratch/rkersten_root/rkersten0/rkersten/Aglaonema_red/entap_config.ini
```
- An entap_config.ini file and a entap_run.params file are required in the EnTAP directory to run EnTAP.
- Generate entap_config.ini:
```
nano entap_config.ini
```
with following input and save:
```
data-generate=false
data-type=0,
entap-db-bin=/sw/pkgs/arc/EnTAP/2.0.0/db/bin/entap_database.bin
entap-graph=/sw/pkgs/arc/EnTAP/2.0.0/bin/entap_graphing.py
rsem-calculate-expression=/sw/pkgs/med/RSEM/1.3.3/bin/rsem-calculate-expression
rsem-sam-validator=/sw/pkgs/med/RSEM/1.3.3/bin/rsem-sam-validator
rsem-prepare-reference=/sw/pkgs/med/RSEM/1.3.3/bin/rsem-prepare-reference
convert-sam-for-rsem=/sw/pkgs/med/RSEM/1.3.3/bin/convert-sam-for-rsem
transdecoder-long-exe=/sw/pkgs/arc/EnTAP/2.0.0/bin/libs/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs
transdecoder-predict-exe=/sw/pkgs/arc/EnTAP/2.0.0/bin/libs/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict
interproscan-exe=interproscan.sh
diamond-exe=/sw/pkgs/arc/EnTAP/1.0.1/bin/libs/bin/diamond
```
- Generate	entap_run.params:
```
fpkm=0.5
align=
transdecoder-m=100
transdecoder-no-refine-starts=false
out-dir=entap_outfiles
overwrite=false
resume=false
input=/path/to/directory/Aglaonema_pearson07.fasta
database=/path/to/database/Uniprot/bin/uniprot_sprot.dmnd
no-trim=false
threads=1
output-format=1,3,4,7,
hgt-donor=
hgt-recipient=
hgt-gff=
ontology_source=0,
eggnog-contaminant=false
eggnog-dbmem=false
eggnog-sensitivity=more-sensitive
interproscan-db=
qcoverage=50
tcoverage=50
e-value=1e-05
uninformative=conserved,predicted,unknown,unnamed,hypothetical,putative,unidentified,uncharacterized,uncultured,uninformative,
diamond-sensitivity=very-sensitive
```
- In the ‘entap_outfiles/final_results’ directory, the entap_results.tsv file was opened in Excel and sorted (A-->Z) to be combined with the correspondingly sorted list of transcripts with Pearson coefficients of 1-0.7 so that annotations were combined with their transcript names and Pearson values.
- The resulting table was searched for target gene families such as 2-oxoglutarate-dependent dioxygenases involved in specialized metabolism.
