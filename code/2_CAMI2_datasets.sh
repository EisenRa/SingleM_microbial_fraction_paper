
### Download a samples from each CAMI2 simulated metagenome type:

cd 1_References/CAMI2

## Marine associated:
# BIOM profiles of a deep-sea environment, using 155 newly sequenced marine isolate genomes from this environment and 622 genomes with matching taxonomic provenance from MarRef
# Additionally, 200 newly sequenced circular elements including plasmids and viruses were added
wget https://frl.publisso.de/data/frl:6425521/marine/short_read/marmgCAMI2_setup.tar.gz
tar -xvzf marmgCAMI2_setup.tar.gz
mv simulation_short_read/ marine
wget https://frl.publisso.de/data/frl:6425521/marine/marmgCAMI2_genomes.tar.gz
tar -xzvf marmgCAMI2_genomes.tar.gz
mv simulation_short_read/genomes/ marine

## Plant associated (90% bacterial, 9% fungi, 1% arabidopsis)
# 894 genomes. Of these, 224 are from the proGenomes terrestrial representative genomes, 216 are newly sequenced genomes from an A. thaliana root rhizosphere, 55 are fungal genomes associated with the rhizosphere, 398 are plasmids or circular elements and one A. thaliana genome
wget https://frl.publisso.de/data/frl:6425521/plant_associated/short_read/rhimgCAMI2_setup.tar.gz
tar -xvzf rhimgCAMI2_setup.tar.gz
mv simulation_short_read plant
wget https://frl.publisso.de/data/frl:6425521/plant_associated/rhimgCAMI2_genomes.tar.gz
tar -xvzf rhimgCAMI2_genomes.tar.gz
mv simulation_short_read/genomes plant

## Strain madness:
# 408 newly sequenced genomes, of which 97% (395) had a closely related strain.
wget https://frl.publisso.de/data/frl:6425521/strain/short_read/strmgCAMI2_setup.tar.gz
tar -xvzf strmgCAMI2_setup.tar.gz
mv short_read strain
wget https://frl.publisso.de/data/frl:6425521/strain/strmgCAMI2_genomes.tar.gz
tar -xvzf strmgCAMI2_genomes.tar.gz
mv short_read/source_genomes/ strain/

# Clean up
rmdir short_read && rmdir simulation_short_read && rm *.gz

## Create coverage profiles to feed into gen-paired-reads
# Add general header information for gen-paired-reads
echo -e '[general]\noutput_sample_name = 2_Simulated_reads/CAMI2_marine0\ninsert_size = 30\ninsert_size_std = 1\nshort_read_length = 150\nerror_rate = 0.005\n' > marine_header.txt
echo -e '[general]\noutput_sample_name = 2_Simulated_reads/CAMI2_plant0\ninsert_size = 30\ninsert_size_std = 1\nshort_read_length = 150\nerror_rate = 0.005\n' > plant_header.txt
echo -e '[general]\noutput_sample_name = 2_Simulated_reads/CAMI2_strain0\ninsert_size = 30\ninsert_size_std = 1\nshort_read_length = 150\nerror_rate = 0.005\n' > strain_header.txt


# Loop over coverage files to get it compatible with gen-paired-reads
while read genome coverage;
  do echo "[$PWD/marine/genomes/$genome]" && echo "coverage = $coverage";
done < marine/coverage_new0.tsv >> marine0_coverage.txt

while read genome coverage;
  do echo "[$PWD/plant/genomes/$genome]" && echo "coverage = $coverage";
done < plant/coverage_new0.tsv >> plant0_coverage.txt

while read genome coverage;
  do echo "[$PWD/strain/source_genomes/$genome]" && echo "coverage = $coverage";
done < strain/coverage_new0.tsv >> strain0_coverage.txt

cat marine_header.txt marine0_coverage.txt > CAMI2_marine0.ini
cat plant_header.txt plant0_coverage.txt > CAMI2_plant0.ini
cat strain_header.txt strain0_coverage.txt > CAMI2_strain0.ini

rm *.txt

##Some refenences have contigs <390 bp (therefore can't use them), fix:

for i in source_genomes*.fasta;
  do reformat.sh \
  in=$i \
  out=${i/.fasta/_FIXED.fasta} \
  minlength=1000;
done

for i in source_genomes/*FIXED.fasta;
  do mv $i ${i/_FIXED/};
done


##Create the simulated metagenomes
conda activate /home/projects/ku-cbd/people/rapeis/.conda-anvio-7.1

#n.b. Neither Phoma_radicina nor Aspergillus_fumigatus exists in the CAMI2 plant
# genomes, despite being in the coverage files. Had to remove it.




for i in *-R1.fastq.gz;
do ../../../../singlem_versions/singlem_dev_03_08_22/bin/singlem pipe \
--singlem_metapackage ../../../0_Database/S3.metapackage_20220513.smpkg/ \
--forward $i \
--reverse ${i/-R1.fastq/-R2.fastq} \
--threads 40 \
--assignment-method naive_then_diamond \
--assignment-singlem-db ../../../0_Database/S3.metapackage_20220513.smpkg/gtdb_r207.reassigned.v5.sdb \
--archive-otu-table ${i/-R1.fastq.gz/_AUGUST_pipe.json} \
--otu-table ${i/-R1.fastq.gz/_AUGUST_pipe.tsv};
done


for i in *.json;
do ../../../../singlem_versions/singlem_dev_03_08_22/bin/singlem condense \
--singlem_metapackage ../../../0_Database/S3.metapackage_20220513.smpkg/ \
--input-archive-otu-table $i \
--output-otu-table ${i/.json/_condense.tsv} \
--apply-expectation-maximisation \
--trim-percent 10;
done



#Calculate mbp of plasmids in plant/marine metagenomes:
filterbyname.sh in=CAMI2_marine0-R1.fastq.gz out=plasmids_marine.fastq names=RNODE include=t substring=t
filterbyname.sh in=CAMI2_plant0-R1 out=plasmids_plant.fastq names=RNODE include=t substring=t









# Reads are interleavened, so split them using bbmap's reformat.sh
for i in *.fastq.gz; do reformat.sh in=$i out1=${i/.fastq/_R1.fastq} out2=${i/.fastq/_R2.fastq}; done

# Run SingleM pipe!
for i in *R1.fastq.gz; do
  singlem pipe \
          --singlem_metapackage ../../2022_SingleM_metagenome/0_Database/S3.0.1.metapackage_20211101.smpkg/ \
          --forward $i \
          --reverse ${i/R1.fastq/R2.fastq} \
          --threads 40 \
          --otu-table $(basename ${i/R1.fastq.gz/pipe.tsv}) \
          --output-extras;
    done

# Run SingleM condense!
for i in *pipe.tsv; do
  singlem condense \
          --singlem_metapackage ../../2022_SingleM_metagenome/0_Database/S3.0.1.metapackage_20211101.smpkg/ \
          --input-otu-tables $i \
          --output-otu-table ${i/_pipe.tsv/_condense.tsv} \
          --trim-percent 10;
    done
