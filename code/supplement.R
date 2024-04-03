mkdir singlem_supplement
cp 3_Outputs/7_Dereplication/animal_coassembly/dereplicated_genomes/*_renamed.fa.gz singlem_supplement/
for i in singlem_supplement/*.gz; do mv $i ${i/_renamed/}; done
cut -f1,2 3_Outputs/8_GTDB-tk/animal_coassembly/animal_coassembly_combined_summary.tsv > singlem_supplement/taxonomy.tsv

#need checkm2 output for supplement
checkm2 predict --threads 8 --input singlem_supplement/*.fa.gz --output-directory singlem_supplement/checkm2_hyena --database_path /projects/ehi/data/0_Environments/databases/CheckM2_database/uniref100.KO.1.dmnd 

singlem supplement \
  --no-dereplication \
  --new-genome-fasta-files *.gz \
  --new-taxonomies taxonomy.tsv \
  --checkm2-quality-file quality_report.tsv \
  --threads 8 \
  --output-metapackage hyena_supplement 