conda install -c conda-forge mamba

mamba env create -f WYNenvironment.yml

mamba init

mamba activate WYNT

wget http://ftp.ebi.ac.uk/pub/databases/ip  d/imgt/hla/hla_gen.fasta -O Data/hla_gen.fasta

makeblastdb -in Data/hla_gen.fasta -out Data/HLAdatabase/HLA_db -dbtype nucl -title HLA_db

