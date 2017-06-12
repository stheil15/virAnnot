# VirAnnot

VirAnnot was build to ease the assembly, blast search and taxonomic annotation of metagenomic NGS data. It is used in Virologie team of [UMR1332 BFP](http://www6.bordeaux-aquitaine.inra.fr/bfp) laboratory at INRA.

It is designed to identify viruses in plants but it can be used to assemble and annotate any sequences with the NCBI taxonomy.

NR and NT must be present localy and NCBI taxonomy is loaded in an SQLITE database with a provided script.

## Prerequisite

#### External programs:
* NCBI Blast+ suite
* SQLite
* Mummer
* Cutadapt ()
* ETE tree
* idba ud
* Bowtie2
* Prinseq


#### External databases:
* NCBI [nr, nt](ftp://ftp.ncbi.nlm.nih.gov/blast/db/)
* NCBI [Taxonomy](ftp://ftp.ncbi.nih.gov/pub/taxonomy)
* [PFAM](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/) (rpsblast files, fasta files, and smp files)

#### Perl external libraries:
* Getopt::Long
* File::Basename
* DBI
* Data::Dumper
* Bioperl
* Color::Rgb
* List::Util
* Excel::Writer

#### Perl included libraries:
* Tools::Fasta
* Tools::Fastq
* Tools::Blast
* Tools::Taxonomy
* Logger::Logger

## Install

Add tools and launchers folders to your $PATH.
```bash
export PATH=/path/to/tools:/path/to/launchers:$PATH
```
Add lib folder to your $PERL5LIB.
```bash
export PERL5LIB=/path/to/lib:$PERL5LIB
```

Download and load NCBI taxonomy into a SQLITE database.
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz ; gunzip taxdump.tar.gz; tar -xf taxdump.tar;
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz ; gunzip prot.accession2taxid.gz;
# Optionally you can combine multiple accession2taxid file with a simple cat. Here dead_prot.accession2taxid and prot.accession2taxid.
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_prot.accession2taxid.gz ; gunzip dead_prot.accession2taxid.gz;
cat prot.accession2taxid dead_prot.accession2taxid > acc2taxid.prot
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz ; gunzip nucl_gb.accession2taxid.gz;
# Optional
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz ; gunzip nucl_wgs.accession2taxid.gz;
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_wgs.accession2taxid.gz ; gunzip dead_wgs.accession2taxid.gz
cat nucl_wgs.accession2taxid nucl_gb.accession2taxid dead_wgs.accession2taxid > acc2taxid.nucl
loadTaxonomy.pl -struct taxonomyStructure.sql -index taxonomyIndex.sql -acc_prot acc2taxid.prot -acc_nucl acc2taxid.nucl -names names.dmp -nodes nodes.dmp
```

Download NCBI nr et nt Blast files.
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz
```
Download PFAM files for RPSBLAST.
```bash
wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Pfam_LE.tar.gz
wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/fasta.tar.gz
wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz
```
Create Taxonomy statistic for each PFAM profile and load it inta the sqlite database.
```bash
\ls -1 *.FASTA | sed 's,^\(.*\)\.FASTA,gi2taxonomy.pl -i & -o \1.tax.txt -r,' | bash
listPath.pl -d . | grep 'tax.txt' > idx
taxo_profile_to_sql.pl -i idx -o taxo_profile.sql
sqlite3 taxonomy.tmp.sqlite < taxo_profile.sql
```
