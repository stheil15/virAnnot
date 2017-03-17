Welcome to VirAnnot's documentation!
====================================

VirAnnot was build to ease the assembly, blast search and taxonomic annotation of metagenomic NGS data. It is used in Virologie team of `UMR1332 BFP <http://www6.bordeaux-aquitaine.inra.fr/bfp>`_ laboratory at INRA.

It was designed to identify viruses in plants but it can be used to assemble and then annotate any sequences with the NCBI taxonomy.

NR and NT must be present localy and NCBI taxonomy is loaded in SQLITE database with a provided script.

.. toctree::
   :maxdepth: 2



Prerequisite
============
External programs
-----------------
* NCBI Blast+ suite
* SQLite
* Mummer
* Blat
* Bowtie2
* Cutadapt ()
* ETE tree

External databases
------------------
* NCBI nr, nt (ftp://ftp.ncbi.nlm.nih.gov/blast/db/)
* NCBI Taxonomy (ftp://ftp.ncbi.nih.gov/pub/taxonomy)
* PFAM (ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/) (rpsblast files, fasta files, and smp files)

Perl external libraries
-----------------------
* Getopt::Long
* File::Basename
* DBI
* Data::Dumper
* Bioperl
* Color::Rgb
* List::Util
* Excel::Writer

Perl included libraries
-----------------------
* Tools::Fasta
* Tools::Fastq
* Tools::Blast
* Tools::Taxonomy
* Logger::Logger

Install
=======
Add tools and launchers folders to your $PATH.







* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
