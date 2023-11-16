## Code for creating simulated Zymo mock communities
export THREADS=40

### Download Zymo mock communities reference genomes
wget https://s3.amazonaws.com/zymo-files/BioPool/ZymoBIOMICS.STD.refseq.v2.zip
unzip ZymoBIOMICS.STD.refseq.v2.zip

#N.B. the Cryptococcus_neoformans_draft_genome is REALLY bad. A lot of contigs are < 1000 bp,
#which means we can't simulate reads from them. I'll filter contigs < 1 kb.
reformat.sh \
in=1_References/ZymoBIOMICS.STD.refseq.v2/Genomes/Cryptococcus_neoformans_draft_genome.fasta \
out=1_References/ZymoBIOMICS.STD.refseq.v2/Genomes/Cryptococcus_neoformans_draft_genome_1kb.fasta \
minlength=1000

rm 1_References/ZymoBIOMICS.STD.refseq.v2/Genomes/Cryptococcus_neoformans_draft_genome.fasta

#Also, headers are not correct for coverM input, so change them for the yeast assemblies:
rename.sh \
in=1_References/ZymoBIOMICS.STD.refseq.v2/Genomes/Cryptococcus_neoformans_draft_genome_1kb.fasta \
out=1_References/ZymoBIOMICS.STD.refseq.v2/Genomes/Cryptococcus_neoformans_draft_genome_1kb_RF.fasta \
prefix=Cryptococcus_

rename.sh \
in=1_References/ZymoBIOMICS.STD.refseq.v2/Genomes/Saccharomyces_cerevisiae_draft_genome.fasta \
out=1_References/ZymoBIOMICS.STD.refseq.v2/Genomes/Saccharomyces_cerevisiae_draft_genome_RF.fasta \
prefix=Saccharomyces_

rm 1_References/ZymoBIOMICS.STD.refseq.v2/Genomes/Saccharomyces_cerevisiae_draft_genome.fasta
rm 1_References/ZymoBIOMICS.STD.refseq.v2/Genomes/Cryptococcus_neoformans_draft_genome_1kb.fasta

#Finally, Bacillus assembly has a funky header, so rename: ">BS.pilon.polished.v3.ST170922"
sed -i'' 's/BS.pilon.polished.v3.ST170922/Bacillus_subtilis_complete_genome/' 1_References/ZymoBIOMICS.STD.refseq.v2/Genomes/Bacillus_subtilis_complete_genome.fasta

### Create simulated reads from these genomes
#https://github.com/merenlab/reads-for-assembly
#N.B., takes uncompressed fastas as input
gen-paired-end-reads 00_Code/Simulated_Zymo.ini
mv 2_Simulated_reads/Simulated_Zymo-R1.fastq 2_Simulated_reads/Simulated_Zymo_R1.fastq
mv 2_Simulated_reads/Simulated_Zymo-R2.fastq 2_Simulated_reads/Simulated_Zymo_R2.fastq

pigz -p $THREADS 2_Simulated_reads/Simulated_Zymo*.fastq

### Spike in simulated human reads into simulated Zymo mock.
#Note, the sim. Zymo mock has 0.38 Gbp, whereas simulated human reads are ~6.6 Gbp
#Need to randomly sample the human reads to 1.14 Gbp -> 25% mock vs 75% human
seqtk sample -s 1337 2_Simulated_reads/Human_reads-R1.fastq.gz 3800000 > 2_Simulated_reads/Human_reads_3-8M_1.fastq
seqtk sample -s 1337 2_Simulated_reads/Human_reads-R2.fastq.gz 3800000 > 2_Simulated_reads/Human_reads_3-8M_2.fastq
pigz -p $THREADS 2_Simulated_reads/Human_reads_3-8M*.fastq

cat 2_Simulated_reads/Human_reads_3-8M_1.fastq.gz 2_Simulated_reads/Simulated_Zymo_R1.fastq.gz > 2_Simulated_reads/Simulated_Zymo_Spiked_Homo_R1.fastq.gz
cat 2_Simulated_reads/Human_reads_3-8M_2.fastq.gz 2_Simulated_reads/Simulated_Zymo_R2.fastq.gz > 2_Simulated_reads/Simulated_Zymo_Spiked_Homo_R2.fastq.gz

#Create two more metagenomes with spiked in arabidopsis and plasmodium DNA
cat 2_Simulated_reads/Arabadopsis_reads-R1.fastq.gz 2_Simulated_reads/Simulated_Zymo_R1.fastq.gz > 2_Simulated_reads/Simulated_Zymo_Spiked_Arabidopsis_R1.fastq.gz
cat 2_Simulated_reads/Arabadopsis_reads-R2.fastq.gz 2_Simulated_reads/Simulated_Zymo_R2.fastq.gz > 2_Simulated_reads/Simulated_Zymo_Spiked_Arabidopsis_R2.fastq.gz
cat 2_Simulated_reads/Plasmodium_reads-R1.fastq.gz 2_Simulated_reads/Simulated_Zymo_R1.fastq.gz > 2_Simulated_reads/Simulated_Zymo_Spiked_Plasmodium_R1.fastq.gz
cat 2_Simulated_reads/Plasmodium_reads-R2.fastq.gz 2_Simulated_reads/Simulated_Zymo_R2.fastq.gz > 2_Simulated_reads/Simulated_Zymo_Spiked_Plasmodium_R2.fastq.gz

### Map simulated reads to reference genomes using Bowtie2
cat 1_References/ZymoBIOMICS.STD.refseq.v2/Genomes/*.fasta > 1_References/ZymoCatted.fna
pigz -p $THREADS 1_References/ZymoCatted.fna

bowtie2-build \
--threads $THREADS \
1_References/ZymoCatted.fna.gz 1_References/ZymoCatted.fna.gz

for i in 2_Simulated_reads/Simulated_Zymo*_R1.fastq.gz; do
  bowtie2 \
          --threads $THREADS \
          -x 1_References/ZymoCatted.fna.gz \
          -1 $i \
          -2 ${i/_R1.fastq/_R2.fastq} \
          --seed 1337 \
          | samtools sort -@ $THREADS -o ${i/_R1.fastq.gz/_Bt2.bam} -;
    done

coverm genome \
        -b 2_Simulated_reads/*.bam \
        -s _ \
        -m relative_abundance \
        -t $THREADS \
        --min-covered-fraction 0 \
        > 3_Outputs/Simulated_zymo_Bt2_coverM.tsv

### Run SingleM pipe on the simulated reads
for i in 2_Simulated_reads/*Zymo*_R1.fastq.gz; do
  singlem pipe \
          --singlem_metapackage 0_Database/S3.0.1.metapackage_20211101.smpkg/ \
          --forward $i \
          --reverse ${i/_R1.fastq/_R2.fastq} \
          --threads $THREADS \
          --include-inserts \
          --otu-table 3_Outputs/1_Simulated_zymo/$(basename ${i/_R1.fastq.gz/_pipe.tsv}) \
          --output-extras;
    done

### Run SingleM condense
for i in 3_Outputs/1_Simulated_zymo/*_pipe.tsv; do
  singlem condense \
          --singlem_metapackage 0_Database/S3.0.1.metapackage_20211101.smpkg/ \
          --input-otu-tables $i \
          --output-otu-table ${i/.tsv/_condense.tsv} \
          --trim-percent 10;
    done


### Updated (Aug) SingleM taxonomy classification method (and databases)
for i in 2_Simulated_reads/*Zymo*_R1.fastq.gz; do
  ../singlem_versions/singlem_dev_03_08_22/bin/singlem pipe \
          --singlem_metapackage 0_Database/S3.metapackage_20220513.smpkg/ \
          --forward $i \
          --reverse ${i/_R1.fastq/_R2.fastq} \
          --threads $THREADS \
          --assignment-singlem-db 0_Database/S3.metapackage_20220513.smpkg/gtdb_r207.reassigned.v5.sdb \
          --assignment-method naive_then_diamond \
          --archive-otu-table 3_Outputs/1_Simulated_zymo/$(basename ${i/_R1.fastq.gz/_AUGUST_pipe.json}) \
          --otu-table 3_Outputs/1_Simulated_zymo/$(basename ${i/_R1.fastq.gz/_AUGUST_pipe.tsv});
    done


for i in 3_Outputs/1_Simulated_zymo/*_AUGUST_pipe.json; do
  ../singlem_versions/singlem_dev_03_08_22/bin/singlem condense \
          --singlem_metapackage 0_Database/S3.metapackage_20220513.smpkg/ \
          --input-archive-otu-table $i \
          --output-otu-table ${i/.json/_condense.tsv} \
          --apply-expectation-maximisation \
          --trim-percent 10;
    done
