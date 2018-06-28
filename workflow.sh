#
# download the genomes
#
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

#
# check the input data
#
<data-tilapia/tilapia_OR_gusta_vomer_TAARs_reference.fasta grep '^>' | sort | uniq -d
# there is one duplicate sequence..

#
# try exonerate with coding2genome model, which should do the
# online translation and comparison in protein space
#
# .. not working as it's usual with biologist tools ;)
# (exhausting 12 GB ram while not searching..)
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

# TODO: merge regions close enough

#
# test gmap for the closer species
#
gmap_build -d r-esox-13 -k 13 -D data-genomes/r-esox data-genomes/r-esox/GCA_000150935.1_ASM15093v1_genomic.fna
gmap -d r-esox-13 -D data-genomes/r-esox -f gff3_gene data-tilapia/tilapia_OR_gusta_vomer_TAARs_reference.fasta > data-results/r-esox.gff3

# simple mapping statistics
<data-results/r-esox.gff3 fgrep "mrna1.exon" | egrep -o "Name=[^; ]*" | sort  -u | wc -l
<data-results/r-esox.gff3 grep -v "^#" | cut -f1 | sort -u | wc -l
<data-results/r-esox.gff3 grep -c "^###"
# 69 query genes hit 18 different scaffolds in 94 gene models

# number of hits in different scaffolds
<data-results/r-esox.gff3 cut -f1 | uniq | grep -v "^#" | sort | uniq -c | sort -rn

#
# the simple way to get spliced coding sequences
#
GENOME=data-genomes/r-esox/GCA_000150935.1_ASM15093v1_genomic.fna

# get all coding sequence (cds) matches
# sort it
# merge overlapping regions, respect strand
# add two more columns to get the strand to the correct column
# get the sequence (reverse complement - strand)
# cluster 'exons' closer than 1000 bases together
# dump each cluster as separate fasta sequence
# wrap the lines at 120 chars
#
<data-results/r-esox.gff3 awk '($3 == "CDS")' |
   sort -k1,1 -k4n,4 |
   bedtools merge -s -i - |
   sed -r 's/([+-])/\t\t\1/' |
   bedtools getfasta -s -bedOut -fi $GENOME -bed - |
   bedtools cluster -d 1000 -s |
   awk -F$'\t' '{
    if ($8 != cur_grp) {
      # finish previous, if any
      if (length(seq)) printf("%d\n%s\n", max, seq);

      # start new
      printf(">%s_%d_", $1, $2);
      cur_grp = $8;
      seq = $7;
      max = $3;
    } else {
      seq = seq $7;
      max = $3;
    }
  }
  END {
    printf("%d\n%s\n", max, seq);
  }' |
  fold -w120 \
> data-results/r-esox.gmap-cds.fna


# not used, replaced by simple awk
#
# now we want to extract exons/CDS of particular matches
# and create special names for the fasta sequences, it's
# easier to use a custom script
#
# conda create -n olfactory
# . activate olfactory
# conda config --add channels bioconda
# conda install pysam


