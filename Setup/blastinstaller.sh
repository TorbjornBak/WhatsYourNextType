#Creating the environment and blast database needed for the pipeline

conda env create -f Setup/WYNenvironment.yml

conda activate WYNT

wget https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_gen.fasta -O Data/hla_gen.fasta

makeblastdb -in Data/hla_gen.fasta -out Data/HLAdatabase/HLA_db -dbtype nucl -title HLA_db
