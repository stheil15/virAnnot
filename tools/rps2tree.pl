#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use Pod::Usage;
use Logger::Logger;
use Tools::Taxonomy;
use Tools::Blast;
use Tools::Fasta;
use DBI;
use Data::Dumper;
use Bio::Seq ;
use Cwd 'abs_path', 'cwd' ;
use Color::Rgb;
use List::Util 'shuffle';
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use GD::Simple;
use Expect;

my $VERSION = '1.3' ;
my $lastmodif = '2018-06-15' ;

my @input_files;
my @id_samples;
my $db = '/media/data/db/taxonomy/taxonomy.sqlite';
my $cdd_fasta_path = '/media/data/db/ncbi/cdd/fasta/';
my $blast_db_path = '/media/data/db/ncbi/nr/nr';
my $blast_type = 'blastx';
my $verbosity=1;
my @seq_fasta;
my @ecsv_files;
my $min_prot_length = 100;
my $viral_portion=0.30;
my $outdir='rps2tree';
my $help;
my $group_file='';
my $p_identity=90;


GetOptions(
  "i|input:s"      => \@input_files,
  "id:s"           => \@id_samples,
  "e|ecsv:s"       => \@ecsv_files,
  "o|outdir:s"     => \$outdir,
  "v|verbosity=i"  => \$verbosity,
  "blast_db=s"     => \$blast_db_path,
  "blast_type=s"   => \$blast_type,
  "db=s"           => \$db,
  "s=s"            => \@seq_fasta,
  "g|group=s"      => \$group_file,
  "mp|min_prot=i"      => \$min_prot_length,
  'vp|viral_portion=f' => \$viral_portion,
  'perc=s'         => \$p_identity,
  "h|help"         => \$help,
);


if($verbosity > 1){
  Logger::Logger->changeMode( $verbosity );
}


&main;


sub main {
  my $self ={};
  bless $self;
  _set_options($self);
  my $cwd = cwd();
  my $rps2tree_folder = $cwd . '/' . $outdir ;
  if(! -e $rps2tree_folder){
    `mkdir $rps2tree_folder`;
  }
  $self->{_colors_array} = _read_rgb_file($self,'/home/stheil/git-repo/virAnnot/lib/rgb.txt');

  $self->{_color_obj} = new Color::Rgb(rgb_txt=>'/home/stheil/git-repo/virAnnot/lib/rgb.txt');
  if($group_file ne ''){
    $logger->info('Treating group file...');
    $self->{_group} = &read_group_file($self,$group_file);
  }
  else{
    for(my $i=0;$i<=$#{$self->{id_samples}};$i++){
      $self->{_group}->{cat_to_color}->{$self->{id_samples}->[$i]} = $self->{_color_obj}->hex($self->{_colors_array}->[$i],'');
      $self->{_group}->{sample_to_color}->{$self->{id_samples}->[$i]} = $self->{_color_obj}->hex($self->{_colors_array}->[$i],'');
    }
  }

  _draw_legend($self);
  _cut_sequences($self);
  _align_with_ref($self);
  $logger->info('Get stats.');
  my $data = _get_global_stats($self,);
  $logger->info('Print excel.');
  _print_excel($self,$outdir . '/rps2tree_stats.xlsx',$data);
  $logger->info('Print csv.');
  _print_csv($self,'cluster_nb_reads.csv',$data);
  $logger->info('Map file.');
  _print_map_file($self,$outdir.'/map.txt',$data);
  _create_html($self,$outdir.'/map.txt',$outdir);
}


sub _draw_legend {
  my ($self)=@_;
  $logger->info('Drawing legend.');
  my $img = GD::Simple->new(400,(scalar(keys(%{$self->{_group}->{cat_to_color}}))*30));
  my $y=10;
  foreach my $grp (keys(%{$self->{_group}->{cat_to_color}})){
    $img->bgcolor($self->{_color_obj}->hex2rgb($self->{_group}->{cat_to_color}->{$grp}));
    $img->rectangle(5,(5+$y),20,(20+$y));
    $img->moveTo(30,(20+$y));
    $img->string($grp);
    $y+=20;
  }
  my $png_file = $outdir . '/legend.png';
  open(PNG,">$png_file") || $logger->logdie('Cannot open file ' . $png_file . ' for writing');
  print PNG $img->png();
  close PNG;
}


sub read_group_file {
  my ($self,$file)=@_;
  my $h;
  open(FILE,$file) || $logger->logdie('Group file not found. '. $file);
  my $j=0;
  while(<FILE>){
    chomp;
    my @l = split(',',$_);
    $h->{cat_to_color}->{$l[0]} = $self->{_color_obj}->hex($self->{_colors_array}->[$j],'');
    # $h->{cat_to_color}->{$l[0]} = $color_names[$j];
    for(my $i=1;$i<=$#l;$i++){
      $h->{sample_to_color}->{$l[$i]} = $self->{_color_obj}->hex($self->{_colors_array}->[$j],'');
      # $h->{sample_to_color}->{$l[$i]} = $color_names[$j];
    }
    $j++;
  }
  close FILE;
  return $h;
}


sub _create_html {
  my ($self,$map_file,$outdir)=@_;
  my $cmd = 'rps2tree_html.py -m ' . $map_file . ' -o ' . $outdir;
  $logger->info($cmd);
  $logger->debug($cmd);
  `$cmd`;
}


sub _print_map_file {
  my ($self,$map_file,$data)=@_;
  open(MAP,">$map_file");
  my @headers = ('cdd_id','align_files','tree_files','cluster_files','cluster_nb_reads_files','pairwise_files','description','full_description');
  print MAP '#' . join("\t",@headers) . "\n";
  foreach my $cdd_id (sort(keys(%{$data}))){
    print MAP $cdd_id ;
    if(defined($self->{_align_files}->{$cdd_id})){
      my @s = split('/',$self->{_align_files}->{$cdd_id});
      print MAP "\t" . $s[$#s-1] . '/' . $s[$#s];
    }
    else{
      print MAP "\t" . '.';
    }

    if(defined($self->{_tree_files}->{$cdd_id})){
      my @s = split('/',$self->{_tree_files}->{$cdd_id});
      print MAP "\t" . $s[$#s-1] . '/' . $s[$#s];
    }
    else{
      print MAP "\t" . '.';
    }

    if(defined($self->{_cluster_files}->{$cdd_id})){
      my @s = split('/',$self->{_cluster_files}->{$cdd_id});
      print MAP "\t" . $s[$#s-1] . '/' . $s[$#s];
    }
    else{
      print MAP "\t" . '.';
    }

    if(defined($self->{_cluster_nb_reads_files}->{$cdd_id})){
      my @s = split('/',$self->{_cluster_nb_reads_files}->{$cdd_id});
      print MAP "\t" . $s[$#s-1] . '/' . $s[$#s];
    }
    else{
      print MAP "\t" . '.';
    }

    if(defined($self->{_pairwise_files}->{$cdd_id})){
      my @s = split('/',$self->{_pairwise_files}->{$cdd_id});
      print MAP "\t" . $s[$#s-1] . '/' . $s[$#s];
    }
    else{
      print MAP "\t" . '.';
    }

    if(defined($self->{cdd_info}->{$cdd_id}->{description})){
      print MAP "\t" . $self->{cdd_info}->{$cdd_id}->{description};
    }
    else{
      print MAP "\t" . '.';
    }

    if(defined($self->{cdd_info}->{$cdd_id}->{full_description})){
      print MAP "\t" . $self->{cdd_info}->{$cdd_id}->{full_description};
    }
    else{
      print MAP "\t" . '.';
    }

    print MAP "\n";
  }
  close MAP;
}


sub _print_excel {
  my ($self,$excel_file,$data)=@_;
  my $workbook = Excel::Writer::XLSX->new( $excel_file );
  foreach my $cdd_id (sort(keys(%{$data}))){
    my $worksheet = $workbook->add_worksheet($cdd_id . '_' . $self->{cdd_info}->{$cdd_id}->{description});
    my @headers = ('#OTU_name');
    push(@headers,sort(keys(%{$data->{$cdd_id}->{sample_present}})));
    $worksheet->write(0,0,\@headers);
    $worksheet->set_column(0,$#headers,12);
    my $row=1;
    foreach my $otu (sort(keys(%{$data->{$cdd_id}->{otu_matrix}}))){
      $worksheet->write($row,0,$otu);
      my $col=1;
      for(my $i=1;$i<=$#headers;$i++){
        if(defined($data->{$cdd_id}->{otu_matrix}->{$otu}->{$headers[$i]})){
          $worksheet->write($row,$col,$data->{$cdd_id}->{otu_matrix}->{$otu}->{$headers[$i]});
        }
        else{
          $worksheet->write($row,$col,'0');
        }
        $col++;
      }
      my @tax = split(';',$data->{$cdd_id}->{otu_annot}->{$otu}->{taxonomy});
      foreach my $t (@tax) {
        $worksheet->write($row,$col,$t);
        $worksheet->set_column($col,$col,length($t));
        $col++;
      }
      $row++;
    }
  }
}


sub _print_csv {
  my ($self,$csv_file,$data)=@_;
  foreach my $cdd_id (sort(keys(%{$data}))){
    my @headers = ('#OTU_name');
    push(@headers,sort(keys(%{$data->{$cdd_id}->{sample_present}})));
    $self->{_cluster_nb_reads_files}->{$cdd_id} = $self->{cdd_folders}->{$cdd_id}.'/'.$csv_file;
    open(CSV,">$self->{cdd_folders}->{$cdd_id}/$csv_file");
    print CSV join("\t",@headers) . "\t" . 'taxonomy' . "\t" . 'seq_list' . "\n";
    foreach my $otu (sort(keys(%{$data->{$cdd_id}->{otu_matrix}}))){
      print CSV $otu;
      for(my $i=1;$i<=$#headers;$i++){
        if(exists($data->{$cdd_id}->{otu_matrix}->{$otu}->{$headers[$i]})){
          print CSV "\t" . $data->{$cdd_id}->{otu_matrix}->{$otu}->{$headers[$i]};
        }
        else{
          print CSV "\t" . '0';
        }
      }
      print CSV "\t" . $data->{$cdd_id}->{otu_annot}->{$otu}->{taxonomy};
      print CSV "\t" . join(',',@{$data->{$cdd_id}->{otu_annot}->{$otu}->{seq_list}});
      print CSV "\n";
    }
    close CSV;
  }
}


sub _get_global_stats {
  my ($self)=@_;
  my $cdd_info;
  foreach my $cdd_id (keys(%{$self->{_cluster_files}})){
    open(CLUST,$self->{_cluster_files}->{$cdd_id});
    while(<CLUST>){
      chomp;
      my @l = split(",",$_);
      for(my $i=2;$i<=$#l;$i++){

        for(my $j=0;$j<=$#{$self->{collection}->{$cdd_id}};$j++){
          # if(defined($self->{collection}->{$cdd_id}->[$j])){
          #   if(defined($self->{collection}->{$cdd_id}->[$j]->{$l[$i]})){
          #     if(defined($self->{_query_annotation}->{$l[$i]})){
          #       if(defined($self->{_query_annotation}->{$l[$i]}->{nb_reads})){
          if(exists($self->{collection}->{$cdd_id}->[$j])){
            if(exists($self->{collection}->{$cdd_id}->[$j]->{$l[$i]})){
              if(exists($self->{_query_annotation}->{$l[$i]})){
                if(exists($self->{_query_annotation}->{$l[$i]}->{nb_reads})){
                  $cdd_info->{$cdd_id}->{otu_matrix}->{$l[0]}->{$self->{id_samples}->[$j]} += $self->{_query_annotation}->{$l[$i]}->{nb_reads};
                }
                else{
                  $cdd_info->{$cdd_id}->{otu_matrix}->{$l[0]}->{$self->{id_samples}->[$j]}++;
                }
                $cdd_info->{$cdd_id}->{sample_present}->{$self->{id_samples}->[$j]}++;
                $cdd_info->{$cdd_id}->{otu_annot}->{$l[0]}->{taxonomy} = $self->{_query_annotation}->{$l[$i]}->{taxonomy};
                push(@{$cdd_info->{$cdd_id}->{otu_annot}->{$l[0]}->{seq_list}},$l[$i]);
                
              }else{
                $cdd_info->{$cdd_id}->{otu_matrix}->{$l[0]}->{$self->{id_samples}->[$i]}++;
                $cdd_info->{$cdd_id}->{otu_annot}->{$l[0]}->{taxonomy} = "not found";
                push(@{$cdd_info->{$cdd_id}->{otu_annot}->{$l[0]}->{seq_list}},$l[$i]);
              }
            }
          }
        }
      }
    }
    close CLUST;
  }
  return $cdd_info;
}


sub _read_rgb_file {
  my ($self,$file)=@_;
  my @a;
  open(RGB,$file);
  while(<RGB>){
    chomp;
    if(/!/){
      next;
    }
    my @line = split("\t",$_);
    push(@a,$line[1]);
  }
  @a = shuffle(@a);
  return \@a;
}


sub _rename_ref_seq {
  my ($self,$cdd_id,$in,$out)=@_;
  if(-e $self->{_cdd_fasta_path} . '/' . $cdd_id . '.FASTA.tax.txt'){
    my $h = {};
    # print $self->{_cdd_fasta_path} . '/' . $cdd_id . '.FASTA.tax.txt' . "\n";
    open(TAX,$self->{_cdd_fasta_path} . '/' . $cdd_id . '.FASTA.tax.txt');
    while(<TAX>){
      chomp;
      my @l = split("\t",$_);
      my @t = split(";",$l[1]);
      $h->{$l[0]} = $t[$#t];
    }
    close TAX;
    open(ORI_FASTA,$in);
    open(CMP_FASTA,">$out");
    while(<ORI_FASTA>){
      chomp;
      if(/^>(\S+)[\s(\S+)]?/){
        my $id = $1;
        if($id =~ /^gi\|(\d+)\|(\S+\|(\S+))\|(.*)$/){
          my $gi = $1;
          my $sp = $2;
          my $acc = $3;
          $acc =~ s/\.\d$//;
          # print CMP_FASTA $_ ;
          if(defined($h->{$acc})){
            $h->{$acc} =~ s/ /_/g;
            $h->{$acc} =~ s/[\(\)]/_/g;
            print CMP_FASTA '>' . $h->{$acc} ;
            if(defined($sp)){
              print CMP_FASTA '_' . $sp;
            }
            print CMP_FASTA "\n";
          }
          else{
            print CMP_FASTA '>' . $sp . "\n";
          }
        }
        else{
          print CMP_FASTA $_ . "\n";
        }
      }
      else{
        print CMP_FASTA $_ . "\n";
      }
      # exit;
    }
    close CMP_FASTA;
    close ORI_FASTA;
  }
}


sub _print_qry_seq {
  my ($self,$cdd_id,$fasta,$config)=@_;
  open(FASTA,">$fasta");
  $self->{_associated_color} = {};
  my $j=0;
  for(my $i =0;$i<=$#{$self->{collection}->{$cdd_id}};$i++){
    # if(defined($self->{collection}->{$cdd_id}->[$i])){
    if(exists($self->{collection}->{$cdd_id}->[$i])){
      foreach my $qry_id (keys(%{$self->{collection}->{$cdd_id}->[$i]})){
        print FASTA '>' . $qry_id . "\n";
        print FASTA $self->{collection}->{$cdd_id}->[$i]->{$qry_id} . "\n";
        $self->{_associated_color}->{$qry_id} = '#' . $self->{_group}->{sample_to_color}->{$self->{id_samples}->[$i]};
      }
      &_print_color_scheme($self,$config);
    }
  }
  close FASTA;
}


sub _seek_ref {
  $logger->info('Seek ref');
  my ($self,$cdd_id,$smp)=@_;
  my $smp_list = $self->{_cdd_folder} . '/smp_list';
  my $ref_pfam = $self->{_cdd_folder} . '/ref_pfam.fasta';
  my $clean_cmd = 'rm ' . $smp_list . "\n";
  open(LIST,">$smp_list");
  print LIST $smp . "\n";
  close LIST;
  my $format_db_cmd = 'makeprofiledb -in ' . $smp_list . ' -out ' . $self->{_cdd_folder} . '/' . $cdd_id . '_rpsdb';
  $clean_cmd .= 'rm ' . $self->{_cdd_folder} . '/' . $cdd_id . '_rpsdb*' . "\n";
  $logger->debug($format_db_cmd);
  `$format_db_cmd`;
  my $ref_acc;
  for(my $i =0;$i<=$#{$self->{collection}->{$cdd_id}};$i++){
    my $blast_hits = $self->{taxoTools}->readCSVextended_regex($self->{ecsvFileList}->[$i],"\t",'Viruses');
    foreach my $q_id (keys(%{$self->{collection}->{$cdd_id}->[$i]})){
      foreach my $m (@{$blast_hits}){
        if($m->{query_id} eq $q_id){
          $self->{_query_annotation}->{$m->{query_id}} = $m;
          $ref_acc->{$m->{accession}}++;
        }
      }
    }
  }
  if(scalar(keys(%{$ref_acc})) > 0){
    my $ref_acc_file = $self->{_cdd_folder} . '/ref_acc.txt';
    $clean_cmd .= 'rm ' . $ref_acc_file . "\n";
    open(REF_ACC,">$ref_acc_file");
    print REF_ACC join("\n",keys(%{$ref_acc}));
    close REF_ACC;
    my $fastacmd = 'fastacmd -d ' . $blast_db_path;
    if($blast_type =~ /^BLASTP$|^BLASTX$/i){
      $fastacmd .= ' -p T ';
    }
    elsif($blast_type =~ /^TBLASTX$|^BLASTN$/i){
      $fastacmd .= ' -p F ';
    }
    else{
    }
    $fastacmd .= ' -i ' . $ref_acc_file . ' > ' . $self->{_cdd_folder} . '/ref.fasta' ;
    $logger->debug($fastacmd);
    `$fastacmd`;

    open(FAS,$self->{_cdd_folder} . '/ref.fasta');
    my $f_line = <FAS>;
    my $lcl=0;
    if($f_line =~ /lcl|/){
      $lcl=1;
    }
    close FAS;

    $clean_cmd .= 'rm ' . $self->{_cdd_folder} . '/ref.fasta' . "\n";
    my $ref_fasta_tool = Tools::Fasta->new('file' => $self->{_cdd_folder} . '/ref.fasta');

    if(scalar(keys(%{$ref_fasta_tool->{index}}))){
      my $rps_cmd = '';
      if($blast_type =~ /^BLASTP$|^BLASTX$/i){
        $rps_cmd = 'rpsblast+';
      }
      elsif($blast_type =~ /^TBLASTX$|^BLASTN$/i){
        $rps_cmd = 'rpstblastn'
      }
      $rps_cmd .= ' -parse_deflines -query ' . $self->{_cdd_folder} . '/ref.fasta' . ' -db ' . $self->{_cdd_folder} . '/' . $cdd_id . '_rpsdb';
      $rps_cmd .= ' -out ' . $self->{_cdd_folder} . '/' . $cdd_id . '_ref.xml -num_threads 5 -max_target_seqs 1 -outfmt 5';
      $logger->debug($rps_cmd);
      `$rps_cmd`;
      $clean_cmd .= 'rm ' . $self->{_cdd_folder} . '/' . $cdd_id . '_ref.xml' . "\n";
      my $rps2ecsv_cmd = 'rps2ecsv.pl -b ' . $self->{_cdd_folder} . '/' . $cdd_id . '_ref.xml' . ' -o ' . $self->{_cdd_folder} . '/' . $cdd_id . '_ref.csv';
      `$rps2ecsv_cmd`;
      $clean_cmd .= 'rm ' . $self->{_cdd_folder} . '/' . $cdd_id . '_ref.csv' . "\n";
      my $ref_pfam_hits = $self->{taxoTools}->readCSVextended($self->{_cdd_folder} . '/' . $cdd_id . '_ref.csv',"\t");
      if(scalar(@{$ref_pfam_hits})){
        open(REF_ACC_CUT,">$ref_pfam");
        for(my $i=0;$i<=$#{$ref_pfam_hits};$i++){
          if(!defined($ref_pfam_hits->[$i]->{startQ}) || ! defined($ref_pfam_hits->[$i]->{endQ})){
            next;
          }
          my $seq='';
          my $acc='';
          my $bioSeq;
          if($lcl){
            $seq = $ref_fasta_tool->retrieveFastaSequence('lcl|'.$ref_pfam_hits->[$i]->{query_id});
            $bioSeq = Bio::Seq->new(-seq => $seq->{'lcl|'.$ref_pfam_hits->[$i]->{query_id}});
            $acc = $ref_pfam_hits->[$i]->{query_id};
          }
          else{
            $seq = $ref_fasta_tool->retrieveFastaSequence($ref_pfam_hits->[$i]->{query_id});
            $bioSeq = Bio::Seq->new(-seq => $seq->{$ref_pfam_hits->[$i]->{query_id}});
            my @a = split('\|',$ref_pfam_hits->[$i]->{query_id});
            $acc = $a[3];
          }
          my $subseq = Bio::Seq->new(-seq => $bioSeq->subseq($ref_pfam_hits->[$i]->{startQ},$ref_pfam_hits->[$i]->{endQ}));
          if($acc =~ /^(\S+)\.\d+$/){
            $acc = $1;
          }
          my $taxid = $self->{taxoTools}->retrieveTaxIdFromAcc($acc,$blast_type);
          my $name = $self->{taxoTools}->retrieveNameFromTaxId($taxid);
          $name =~ s/\s/_/g;
          $name =~ s/[\(\)]/_/g;

          if(length($subseq->seq) >= $self->{_min_prot_length}){
            print REF_ACC_CUT '>' . $name . '_' . $acc . "\n";
            if($blast_type =~ /^BLASTP$|^BLASTX$/i){
              print REF_ACC_CUT $subseq->seq . "\n";
            }
            elsif($blast_type =~ /^TBLASTX$|^BLASTN$/i){
              if($ref_pfam_hits->[$i]->{frame} < 0){
                $subseq = $subseq->revcom;
              }
              my $prot = $subseq->translate();
              print REF_ACC_CUT $prot->seq . "\n";
            }
          }
        }
        close REF_ACC_CUT;
      }
    }
  }
  if($verbosity < 3){
    $logger->debug($clean_cmd);
    `$clean_cmd`;
  }
}


sub _align_with_ref {
  $logger->info('Align with ref :');
  my ($self)=@_;
  my $cwd = cwd();
  my $rps2tree_folder = $cwd . '/' . $outdir ;

  my $clean_cmd;
  foreach my $cdd_id (keys(%{$self->{collection}})){
    $logger->info('---- ' . $cdd_id . ' ----');
    if(-e $self->{_cdd_fasta_path} . '/' . $cdd_id . '.FASTA'){
      my $ori_cdd_fasta = $self->{_cdd_fasta_path} . '/' . $cdd_id . '.FASTA';
      $self->{_cdd_folder} = $rps2tree_folder . '/' . $cdd_id . '_' . $self->{cdd_info}->{$cdd_id}->{description};
      $self->{cdd_folders}->{$cdd_id} = $self->{_cdd_folder};
      if(! -e $self->{_cdd_folder}){
        `mkdir $self->{_cdd_folder}`;
      }
      my $cp_cdd_fasta = $self->{_cdd_folder} . '/' . $cdd_id . '.fasta';
      my $smp_file = $self->{_cdd_fasta_path} . '/' . $cdd_id . '.smp';
      if(! -e $cp_cdd_fasta){
        my $cp_cmd = 'cp ' . $ori_cdd_fasta . ' ' . $cp_cdd_fasta;
        `$cp_cmd`;
      }

      my $cmp_fasta = $self->{_cdd_folder} . '/' . $cdd_id . '.newid.fasta';
      if(! -e $cmp_fasta){
        _rename_ref_seq($self,$cdd_id,$cp_cdd_fasta,$cmp_fasta);
      }
      $clean_cmd .= 'rm ' . $cp_cdd_fasta . "\n";



      my $query_seq_fasta = $self->{_cdd_folder} . '/' . 'qry_seq.fa' ;
      my $config_file = $self->{_cdd_folder} . '/' . 'tree_config.txt';
      if(! -e $query_seq_fasta){
        _print_qry_seq($self,$cdd_id,$query_seq_fasta,$config_file);
      }

      my $seq_to_align = $self->{_cdd_folder} . '/all_seq_to_align.fa';
      _seek_ref($self,$cdd_id,$smp_file);
      # if(! -e $seq_to_align){
      # }

      if(-e $self->{_cdd_folder} . '/ref_pfam.fasta'){
        my $cat_cmd = 'cat ' . $self->{_cdd_folder} . '/ref_pfam.fasta' . ' ' . $query_seq_fasta . ' ' . $cmp_fasta . ' > ' . $seq_to_align;
        `$cat_cmd`;
      }
      else{
        my $cat_cmd = 'cat ' . $query_seq_fasta . ' ' . $cmp_fasta . ' > ' . $seq_to_align;
        `$cat_cmd`;
      }

      my $align_fasta = $self->{_cdd_folder} . '/' . $cdd_id . '_w_qry_aligned.fasta';
      my $tree_file = $self->{_cdd_folder} . '/' . $cdd_id . '_' . $self->{cdd_info}->{$cdd_id}->{description} . '.tree.nwx';
      if(! -e $align_fasta){
        my $ete_cmd = 'ete3 build -w standard_trimmed_fasttree -a ' . $seq_to_align . ' -o ' . $self->{_cdd_folder} . '/ete3 --rename-dup-seqnames' . "\n";
        $logger->debug($ete_cmd);
        $logger->info($ete_cmd);
        `$ete_cmd`;

        my $mv_ete_cmd = 'mv ' . $self->{_cdd_folder} . '/ete3/clustalo_default-trimal01-none-fasttree_full/all_seq_to_align.fa.final_tree.fa' . ' ' . $align_fasta . "\n";
        $mv_ete_cmd .= 'mv ' . $self->{_cdd_folder} . '/ete3/clustalo_default-trimal01-none-fasttree_full/all_seq_to_align.fa.final_tree.nwx' . ' ' . $tree_file . "\n";
        $clean_cmd .= 'rm -r ' . $self->{_cdd_folder} . '/ete3' . "\n";
        $logger->debug($mv_ete_cmd);
        `$mv_ete_cmd`;
        $self->{_align_files}->{$cdd_id} = $align_fasta;
      }
      else{
        $self->{_align_files}->{$cdd_id} = $align_fasta;
      }

      my $matrix_file = $self->{_cdd_folder} . '/' . 'identity_matrix.csv';

      if(! -e $tree_file . '.png'){
        my $ete_tree_cmd = 'ete_tree.py' . ' -f ' . $align_fasta . ' -t ' . $tree_file . ' -c ' . $config_file . ' -m ' . $matrix_file . ' -o ' . $tree_file . '.png';
        $logger->debug($ete_tree_cmd);
        $logger->info($ete_tree_cmd);
        `$ete_tree_cmd`;
        if(-e $tree_file . '.png'){
          $self->{_tree_files}->{$cdd_id} = $tree_file . '.png';
        }
      }
      else{
        $self->{_tree_files}->{$cdd_id} = $tree_file . '.png';
      }
      my $cluster_file = $self->{_cdd_folder} . '/' . 'clusters.csv';
      my $temp_otu_file = $self->{_cdd_folder} . '/' . 'temp_otu_file.csv';
      if(! -e $cluster_file && -e $align_fasta){
        &_compute_pairwise_distance($self,$align_fasta,$matrix_file,$cluster_file,0,$self->{p_identity}, $temp_otu_file);
        $self->{_cluster_files}->{$cdd_id} = $cluster_file;
        $self->{_pairwise_files}->{$cdd_id} = $matrix_file;
      }
      else{
        $self->{_cluster_files}->{$cdd_id} = $cluster_file;
        $self->{_pairwise_files}->{$cdd_id} = $matrix_file;
      }

    }
    $logger->debug($clean_cmd);
  }
  # if($verbosity < 3){
  #   `$clean_cmd`
  # }
}


sub _compute_pairwise_distance {
  # { data , aligned sequences, output cluster file, take gap into account, identity percentage, OTU temporary output file}
  my ($self,$fasta,$matrix_file,$cluster_file,$gap,$hsppercid_opt,$temp_otu_file)=@_;
  $logger->info("Computing pairwise distance with " . $hsppercid_opt . " percentage of identity...");
  open(FASTA,$fasta);
  my $id='';
  my $seq='';
  my $sequences={};
  $logger->info("Reading fasta file...");
  while(<FASTA>){
      chomp;
      if(/>(\S+)/){
          if($seq ne ''){
              my @tab = split(//,$seq);
              $sequences->{$id} = \@tab;
              $seq='';
          }
          $id = $1;
      }
      else{
          $seq .= $_;
      }
  }
  close FASTA;
  my @tab = split(//,$seq);
  $sequences->{$id} = \@tab;
  $logger->debug(scalar(keys(%{$sequences})) . ' loaded.');
  my @matrix;

  my %hash_family ;
  my %array_family ;
  my $kfamily = 0 ;
  my $num_results ;
  my $count = 0;
  my %h = ();

  # CALCULATE THE PERCENT IDENTITY BETWEEN EACH PAIR OF SEQUENCES:
  my @keys1 = sort(keys(%{$sequences})) ; # SORTING THE KEYS MAKES IT EASIER TO TEST
  $logger->info("Calculate identity matrix...");
  for(my $i=0; $i<=$#keys1; $i++){
    my $seq1 = $sequences->{$keys1[$i]};
    my $fam_nb = $kfamily ;
    if ($hash_family{$keys1[$i]}) {
        $fam_nb = $hash_family{$keys1[$i]} ;
        
    }
    else {
            $kfamily += 1 ;
            $fam_nb = $kfamily ;
            $hash_family{$keys1[$i]} = $fam_nb ;
            push (@{$array_family{$fam_nb}}, $keys1[$i]) ;
    }
    for(my $j=0; $j<=$#keys1; $j++){
      my $seq2 = $sequences->{$keys1[$j]};

      my $identic=0;
      my $compared=0;
      for(my $k=0; $k < scalar(@{$seq1}); $k++){
          # mutation, next
          if($seq1->[$k] eq 'X' || $seq2->[$k] eq 'X'){next;}
          # gap in both sequences, next
          if($seq1->[$k] eq '-' && $seq2->[$k] eq '-'){next;}
          # gap in one of the sequence, next
          if($seq1->[$k] eq '-' || $seq2->[$k] eq '-'){next;}
          # identity
          if($seq1->[$k] =~ /$seq2->[$k]/i){
            $identic++;
          }
          $compared++;
      }
      my $percentIdentity;
      # set minimum overlap to 20
      if($compared == 0 || $compared < 20){
        $percentIdentity=0;
      }
      else{
        $percentIdentity = ( $identic / $compared )*100 ;
      }
      $matrix[$i][$j] = $percentIdentity;
      if ($percentIdentity > $hsppercid_opt) {
        $a = [ $keys1[$i], $keys1[$j] ];
        $h{$count} = $a ;
        $count += 1;
      }
      $hash_family{$keys1[$j]} and next ;
      next if ($percentIdentity < $hsppercid_opt );
      $hash_family{$keys1[$j]} = $fam_nb ;
      push (@{$array_family{$fam_nb}}, $keys1[$j]) ;
      
    }
  }
  # Write matrix to file
  open(MAT,">$matrix_file");
  # print MAT "\t" . scalar(@keys1) . "\n" ;
  for(my $i=0; $i<=$#matrix; $i++){
    my $name = $keys1[$i];
    my $dim = scalar(@keys1);
    my $count = 1;
    printf MAT "%50s",$name ;
    # for(my $j=0; $j<=$i; $j++){
    for(my $j=0; $j<=$dim; $j++){
      if($j < $#{$matrix[$i]}){
        if ($count != 0){
          printf MAT ",%.4f", $matrix[$i][$j];
        }else{
          printf MAT ",%.4f", 0;
        }
        # if ($matrix[$i][$j] == 100){
        #   $count = 0;
        # }
      }else{
        # if ( $matrix[$i][$j] == 100) {
          printf MAT ",%.4f", $matrix[$i][$j];
        # }
        # else{
        #   printf MAT ",%.4f", 0;
        # }
      }
    }
    print MAT ",\n";
  }
  close MAT;
  # Run OTUs assignation
  $logger->info("Seek OTUs...");
  `seek_otu.R $matrix_file $cluster_file $hsppercid_opt`;
  print STDERR "# $kfamily families found for " . scalar(@keys1) . " queries\n" ;
}


sub _print_color_scheme {
  my ($self,$config_file)=@_;
  open(CONF,">>$config_file");
  foreach my $qry_id (keys(%{$self->{_associated_color}})){
    print CONF $qry_id . "\t" . $self->{_associated_color}->{$qry_id} . "\n";
  }
  close CONF;
}


sub _cut_sequences {
  $logger->info('Cut sequences.');
  my ($self)=@_;
  for(my $i=0;$i<=$#{$self->{filesList}};$i++){
    my $pfam_hits = $self->{taxoTools}->readCSVextended_regex($self->{filesList}->[$i],"\t",'Viruses');
    my $fasta_tool = Tools::Fasta->new('file' => $self->{seqFileList}->[$i]);
    $logger->info('Csv read');
    # for(my $j=0;$j<=$#{$self->{_pfam_hits}};$j++){
    for(my $j=0;$j<=$#{$pfam_hits};$j++){
      if(!defined($pfam_hits->[$j]->{superkingdom})){next;}
      if($pfam_hits->[$j]->{superkingdom} =~ /Viruses\((\d+\.\d+)\);/){
        my $vir_percent = $1;
        if($vir_percent >= $self->{_viral_portion}){
          my $seq = $fasta_tool->retrieveFastaSequence($pfam_hits->[$j]->{query_id});
          my $seq_length = length($seq->{$pfam_hits->[$j]->{query_id}});
          my $bioSeq = Bio::Seq->new(-seq => $seq->{$pfam_hits->[$j]->{query_id}});
          $logger->info($pfam_hits->[$j]->{startQ} . " " . $pfam_hits->[$j]->{endQ} . " " . $seq_length);
          my $subseq;
          if($pfam_hits->[$j]->{endQ} < $seq_length){
            $subseq = Bio::Seq->new(-seq => $bioSeq->subseq($pfam_hits->[$j]->{startQ},$pfam_hits->[$j]->{endQ}));
          }else{
            $logger->info('seq length');
            $subseq = Bio::Seq->new(-seq => $bioSeq->subseq($pfam_hits->[$j]->{startQ},$seq_length));
          }
          if($pfam_hits->[$j]->{frame} < 0){
            $subseq = $subseq->revcom;
          }
          my $prot = $subseq->translate();
          if(length($prot->seq) > $self->{_min_prot_length}){
            $self->{collection}->{$pfam_hits->[$j]->{cdd_id}}->[$i]->{$pfam_hits->[$j]->{query_id}} = $prot->seq;
            $self->{cdd_info}->{$pfam_hits->[$j]->{cdd_id}}->{full_description} = $pfam_hits->[$j]->{description};
            my @a = split(',',$pfam_hits->[$j]->{description});
            $self->{cdd_info}->{$pfam_hits->[$j]->{cdd_id}}->{description} = $a[1];
            $self->{cdd_info}->{$pfam_hits->[$j]->{cdd_id}}->{description} =~ s/^\s//;
          }
        }
      }
    }
  }
  $logger->info('End Cut sequences.');
}


sub _set_options {
  my ($self)=@_;
  if(scalar(@input_files) > 0){
    if(scalar(@id_samples) == 0){
      foreach my $file (@input_files){
        my($filename, $dirs, $suffix) = fileparse($file);
        my @p = split(/\./,$filename);
        push(@id_samples,$p[0]);
      }
    }
    $self->{id_samples} = \@id_samples;
    # RPS blast files
    $self->{filesList} = \@input_files;
    # contigs
    if(scalar(@seq_fasta) > 0){
      $self->{seqFileList} = \@seq_fasta;
    }
    else{
      $logger->logdie('You must provide sequence files associated with csv files.')
    }
    if(scalar(@ecsv_files) > 0){
      # Blastx files
      $self->{ecsvFileList} = \@ecsv_files;
    }
    else{
      $logger->logdie('You must provide annotation files associated with csv files.')
    }
  }
  else{
    $logger->error('You must provide at least one csv file.');
    &help;
  }
  my %taxonomyParams;
  if(! -e $db){
    $logger->error('Taxonomy sqlite database not found. ' . $db);
    exit;
  }
  else{
    $taxonomyParams{'dbh'} = DBI->connect("dbi:SQLite:dbname=$db","","");
  }
  $self->{taxoTools} = Tools::Taxonomy->new(%taxonomyParams);
  if(-e $cdd_fasta_path){
    $self->{_cdd_fasta_path} = $cdd_fasta_path;
  }
  else{
    $logger->error('Path for CDD fasta not found.');
    exit;
  }
  if($min_prot_length > 0){
    $self->{_min_prot_length} = $min_prot_length;
  }
  else{
    $logger->logdie('Minimum protein length must be > 0.');
  }
  if($viral_portion >= 0 && $viral_portion <= 1){
    $self->{_viral_portion} = $viral_portion;
  }
  if($p_identity > 0 && $viral_portion < 100){
    $self->{p_identity} = $p_identity;
  }
  else{
    $logger->logdie('Viral portion must be set between 0 and 1.');
  }
}


sub help {
my $prog = basename($0) ;
print STDERR <<EOF ;
### $prog $VERSION ###
#
# AUTHOR:     Sebastien THEIL
# VERSION:    $VERSION
# LAST MODIF: $lastmodif
# PURPOSE:    
#

USAGE: perl $prog -i rspblast_csv_1 -id id_sample1 -e blast_cvs_1 ... -i rpsblast_csv_n -id id_sampleN -e blast_csv_n [OPTIONS]

       ### OPTIONS ###
       -i|input           <RPSBLAST CSV>  Rpsblast CSV file from rps2csv.pl script.
       -id                Sample ID.
       -e|ecsv            <BLAST CSV> Blast CSV with taxonomy field from blast2ecsv.pl script.
       -o|outdir          Output directory.
       -s                 Sequences in fasta.
       -mp|min_prot       Minimum query protein length.
       -vp|viral_portion  Minimun portion of viral sequences in PFAM domain to be included.
       -perc              Minimum percentage of identity between two sequences
       -db                NCBI Taxonomy SqliteDB.
       --blast_db         Blast database path.
       --blast_type       Blast program type.
       -g|group           Group file. CSV file, comma separated, one line per group, first column for group name, then samples names. Sequence ID format has to be "SAMPLE_QUERYID".
       -v|verbosity       1 -> 3

       -help|h            Print this help and exit
EOF
exit(1) ;
}
