# download the genomes
# medaka
# is in ensembl
wget ftp://ftp.ensembl.org/pub/release-92/fasta/oryzias_latipes/pep/Oryzias_latipes.MEDAKA1.pep.abinitio.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/oryzias_latipes/pep/Oryzias_latipes.MEDAKA1.pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/oryzias_latipes/dna/Oryzias_latipes.MEDAKA1.dna_sm.toplevel.fa.gz

# R. esox
# the worst genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/150/935/GCA_000150935.1_ASM15093v1/GCA_000150935.1_ASM15093v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/150/935/GCA_000150935.1_ASM15093v1/GCA_000150935.1_ASM15093v1_genomic.gff.gz

# Scleropages
# has refseq build at least
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/624/265/GCF_001624265.1_ASM162426v1/GCF_001624265.1_ASM162426v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/624/265/GCF_001624265.1_ASM162426v1/GCF_001624265.1_ASM162426v1_genomic.gff.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/624/265/GCF_001624265.1_ASM162426v1/GCF_001624265.1_ASM162426v1_protein.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Scleropages_formosus/latest_assembly_versions/GCF_001624265.1_ASM162426v1/GCF_001624265.1_ASM162426v1_translated_cds.faa.gz

# shit...not working as usually
QUERY=data-tilapia/tilapia_OR_gusta_vomer_TAARs_reference.fasta
TARGET=data-genomes/r-esox/GCA_000150935.1_ASM15093v1_genomic.fna
exonerate \
  -M 5000 \
  -Q dna -q $QUERY \
  -T dna -t $TARGET --softmasktarget y \
  -m coding2genome -g y -p blosum62 \
  --percent 90 --showtargetgff y \
> data-results/r-esox.gff

# smaller sample, get timing for 20 sequences, multicore
# <( ) does not work..
<$QUERY head -40 > data-tilapia/q-sample-20.fna
time exonerate \
  -M 5000 --cores 6 \
  -Q dna -q data-tilapia/q-sample-20.fna \
  -T dna -t $TARGET --softmasktarget y \
  -m est2genome -g y -p blosum62 \
  --percent 90 --showtargetgff y -V 50 \
2> data-results/r-esox.out

# this worked a little
# 8 minutes, 1 hit
# vulgar: NC_022205_LG7TAAR_1like__LOC102080301_CDS 0 981 + ABPN01018094.1 341 1304 + 4433 M 663 663 G 3 0 M 132 132 G 15 0 C 3 3 M 165 165
time exonerate   -M 5000 --cores 6  -q data-tilapia/q-sample-20.fna   -t $TARGET --softmasktarget y   -m cdna2genome -g y -p blosum62   --percent 90 --showtargetgff y -V 50 2> data-results/r-esox.out


#
# lastz target query
#
alias lastz='/opt/lastz/bin/lastz'

# simple test
QTEST=data-tilapia/q-sample-20.fna
lastz "$TARGET[multiple]" $QTEST --format=sam- | less -S

# full run
time lastz "$TARGET[multiple]" $QUERY --format=sam > data-results/r-esox.sam
# 16 seconds, 120 sequences matched

time lastz "$TARGET[multiple]" $QUERY --format=general > data-results/r-esox.tsv


