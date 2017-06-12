#!/usr/bin/perl
use strict ;
use File::Basename;
use Data::Dumper;
use Tools::Fasta;
use Tools::Fastq;
use Tools::Taxonomy;
use Logger::Logger;
use Getopt::Long;
use Cwd;
use Cwd 'abs_path', 'cwd';

my $blast_ecsv_file='';
my $selected_taxo_rank=2;
my $seek_up_or_down='up';
my $out_path='assembly_by_taxon';
my $sample_id='';
my $bam_file='';
my $output='abt_scaffolds.fasta';
my $pair_R1='';
my $pair_R2='';
my $singletons='';
my $read_norm=0;
my $mega_merge=0;
my $verbosity=1;
my $n_cpu=10;
my $VERSION = '1.2';
my $lastmodif = '2016-10-27';



GetOptions(
	"ecsv|e=s"         => \$blast_ecsv_file,
	"bam|b=s"          => \$bam_file,
	"s=s"              => \$sample_id,
	"o=s"              => \$output,
	"1=s"              => \$pair_R1,
	"2=s"              => \$pair_R2,
  "rn"               => \$read_norm,
  "seek=s"           => \$seek_up_or_down,
  "singletons=s"     => \$singletons,
  "mm|megamerge"     => \$mega_merge,
  "n_cpu=i"          => \$n_cpu,
	"v|verbosity=i"    => \$verbosity
);

if($verbosity > 1){
    Logger::Logger->changeMode($verbosity);
}


&main;


sub main {
	my $self={};
	bless $self;
	_set_options($self);
	$self->{taxoTools} = Tools::Taxonomy->new();
	$self->{_blast_annotation} = $self->{taxoTools}->readCSVextended($self->{_blast_ecsv_file},"\t");
	_sort_by_rank($self,$self->{_selected_taxo_rank},$self->{_seek_up_or_down});

	_dispatch_reads_from_bam($self);
  _build_assembly_cmd_fastq($self);

  if($self->{_mega_merging}==1){
    _launch_mega_merge($self);
    my $mega_merging_path = cwd() . '/' . $out_path . '/mega_merge/MergedContigs.fasta';
    my @a = split('/',cwd());
    my $cmd = 'sed -i \'s,^>Contig_\([0-9]*\),>' . substr($a[$#a],0,8) . '_\1,\' ' . $mega_merging_path . "\n";
    $cmd .= 'mv ' . $mega_merging_path . ' ' . $self->{_output} ;
    $logger->debug($cmd);
    system($cmd);
  }
  else{
    foreach my $rank (keys(%{$self->{_assembly_info}})){
  		if(defined($self->{_assembly_info}->{$rank}->{scaffolds})){
  			my $cmd = 'cat ' . $self->{_assembly_info}->{$rank}->{scaffolds} . ' >> ' . $self->{_output};
  			$logger->debug($cmd);
  			`$cmd`;
  		}
  	}
  }
}


sub _launch_mega_merge {
  my ($self)=@_;
  my $mega_merging_path = cwd() . '/' . $out_path . '/mega_merge';
  my $cmd='';
  $cmd = 'MeGAMerge-1.1.pl ' . $mega_merging_path ;
  foreach my $rank (keys(%{$self->{_assembly_info}})){
    if(!defined($self->{_assembly_info}->{$rank}->{scaffolds})){
      $logger->warn($rank . ' seems to have failed.');
    }
    else{
      $cmd .= ' ' . $self->{_assembly_info}->{$rank}->{scaffolds};
    }
  }
  $logger->debug($cmd);
  `$cmd`;
}


sub _dispatch_reads_from_bam {
	my ($self)=@_;
	foreach my $rank (keys(%{$self->{_scaffold_by_taxo_rank}})){
    my $rank_path = cwd() . '/' . $out_path . '/' . $rank;

    if(! -e $rank_path){
      `mkdir -p $rank_path`;
    }
		$logger->info($rank . ' containing ' . scalar(@{$self->{_scaffold_by_taxo_rank}->{$rank}}) . ' scaffolds.');
    my $cmd = '';
    if(-e $rank_path . '/read_id.txt'){
      `rm $rank_path/read_id.txt`
    }
    if(scalar(@{$self->{_scaffold_by_taxo_rank}->{$rank}}) > 10000){
      while( my @splice = splice(@{$self->{_scaffold_by_taxo_rank}->{$rank}},0,10000)  ){
        $cmd = 'samtools view -F 4 ' . $self->{_bam_file} . ' ' . join(' ',@splice) . ' | cut -f 1 | sort -u >> ' . $rank_path . '/read_id.txt' . "\n";
        $logger->debug($cmd);
        `$cmd`;
      }
    }
    else{
      $cmd = 'samtools view -F 4 ' . $self->{_bam_file} . ' ' . join(' ',@{$self->{_scaffold_by_taxo_rank}->{$rank}}) . ' | cut -f 1 | sort -u > ' . $rank_path . '/read_id.txt' . "\n";
      $logger->debug($cmd);
      `$cmd`;
    }
    $cmd = 'wc -l ' . $rank_path . '/read_id.txt' . "\n";
    $logger->debug($cmd);
		my $res = `$cmd`;
    chomp $res;
		$logger->debug((split(" ",$res))[0] . ' read ids have been found.');

    $cmd = 'seqtk subseq ' . $self->{_pair_R1} . ' ' . $rank_path . '/read_id.txt' . ' > ' . $rank_path . '/' . $rank . '_R1.fastq' . "\n";
		$cmd .= 'sed -i \'s,^@.*$,& 1:N:0:ACAG,\' ' . $rank_path . '/' . $rank . '_R1.fastq';
    $logger->debug($cmd);
    `$cmd`;

    $cmd = 'seqtk subseq ' . $self->{_pair_R2} . ' ' . $rank_path . '/read_id.txt' . ' > ' . $rank_path . '/' . $rank . '_R2.fastq' . "\n";
		$cmd .= 'sed -i \'s,^@.*$,& 2:N:0:ACAG,\' ' . $rank_path . '/' . $rank . '_R2.fastq';
		$logger->debug($cmd);
    `$cmd`;

    if($singletons ne ''){
      $cmd = 'seqtk subseq ' . $self->{_singletons} . ' ' . $rank_path . '/read_id.txt' . ' > ' . $rank_path . '/' . $rank . '_singletons.fastq' . "\n";
  		$cmd .= 'sed -i \'s,^@.*$,& 2:N:0:ACAG,\' ' . $rank_path . '/' . $rank . '_singletons.fastq';
  		$logger->debug($cmd);
      `$cmd`;
    }
		# my @a = map {$_ . '/1'} @read_list;
		# print PAIR_1 $pair_R1_index->retrieveFastqBlock(\@a);
		# @a = map {$_ . '/2'} @read_list;
		# print PAIR_2 $pair_R2_index->retrieveFastqBlock(\@a);
	}
}


sub _build_assembly_cmd_fastq {
	my ($self)=@_;
  foreach my $rank (keys(%{$self->{_scaffold_by_taxo_rank}})){
    $logger->debug('########## ' . $rank . ' ##########');
  	my $rank_path = cwd() . '/' . $out_path . '/' . $rank;
    my $final_pairs = $rank . '_merge.fastq';
    my $final_single = '';
  	my $cmd = '';
    $cmd .= 'cd ' . $rank_path . "\n";
    $cmd .= 'fastqutils merge ' . $rank . '_R1.fastq' . ' ' . $rank . '_R2.fastq' . ' > ' . $final_pairs . "\n";
    $logger->debug($cmd);
    # `$cmd`;
    if($self->{_read_norm}==1){
      $cmd = 'cd ' . $rank_path . "\n";
      $cmd .= 'normalize-by-median.py -N 12 -k 20 -C 40 -x 5e8 -R ' . 'khmer_norm.log'  . ' -s ' . 'c40k20norm.kh' . ' -p ' . $final_pairs . "\n";
      $cmd .= 'filter-abund.py -V ' . 'c40k20norm.kh' . ' ' . $rank . '_merge.fastq.keep' . "\n";
      $cmd .= 'extract-paired-reads.py ' . ' --output-paired ' . $rank . '_final_pairs.fastq' . ' --output-single ' . $rank . '_final_single.fastq' . ' ' . $rank . '_merge.fastq.keep.abundfilt' . "\n";
      $logger->debug($cmd);
      # `$cmd`;
      $final_pairs = $rank . '_final_pairs.fastq';
      $final_single = $rank . '_final_single.fastq';
      if(! -z $final_pairs){
        $cmd ='cd ' . $rank_path . "\n";
        $cmd .= 'metaspades.py';
        $cmd .= ' --12 ' . $rank . '_final_pairs.fastq';
        if(! -z $final_single){
          $cmd .= ' -s ' . $final_single;
        }
        $cmd .= ' -o ' . 'spades' ;
      	$cmd .= ' -t ' . $self->{_n_cpu};
        $logger->debug($cmd);

      }
      else{
        $logger->warn('No pairs found for ' . $rank);
      }
    }
    else{
      $cmd ='cd ' . $rank_path . "\n";
      $cmd .= 'metaspades.py';
      if(! -z $rank . '_R1.fastq' && ! -z $rank . '_R2.fastq'){
        $cmd .= ' -1 ' . $rank . '_R1.fastq';
        $cmd .= ' -2 ' . $rank . '_R2.fastq';
      }
      if(! -z $rank .'_singletons.fastq'){
        $cmd .= ' -s ' . $rank .'_singletons.fastq';
      }
      $cmd .= ' -o ' . 'spades' ;
      $cmd .= ' -t ' . $self->{_n_cpu};
      $logger->debug($cmd);
    }
    $self->{_cmd}->{$rank} = $cmd;
  }
}


sub _launch_cmd {
	my ($self)=@_;
	foreach my $rank (keys(%{$self->{_cmd}})){
		$logger->info('Assembly for '. $rank);
		my $res = `$self->{_cmd}->{$rank}`;
		if(_parse_spades_exec_log($self,$res)){
			$logger->warn('Assembly rank ' . $rank . ' failed. Retry with shorter kmer..');
			$res = `$self->{_cmd}->{$rank} -k 21,33`;
			if(-e 'core'){
				`rm core`;
			}
			if(_parse_spades_exec_log($self,$res)){
				$logger->warn('Assembly rank ' . $rank . ' failed two times. Next rank..');
				if(-e 'core'){
					`rm core`;
				}
				next;
			}
			$self->{_assembly_info}->{$rank} = _get_assembly_info($self,$res);
		}
		else{
			$self->{_assembly_info}->{$rank} = _get_assembly_info($self,$res);
		}
		if(defined($self->{_assembly_info}->{$rank}->{scaffolds})){
			_rename_scaff($self,$rank,$self->{_assembly_info}->{$rank}->{scaffolds});
		}
	}
	if(-e $self->{_output}){
		`rm $self->{_output}`;
	}
	foreach my $rank (keys(%{$self->{_assembly_info}})){
		if(defined($self->{_assembly_info}->{$rank}->{scaffolds})){
			my $cmd = 'cat ' . $self->{_assembly_info}->{$rank}->{scaffolds} . ' >> ' . $self->{_output};
			$logger->debug($cmd);
			`$cmd`;
		}
	}
}


sub _rename_scaff {
	my ($self,$rank,$scaff_path)=@_;
	my $cmd = 'sed -i \'s,^>NODE_\([0-9]*\)_\(.*\)$,>' . substr($rank,0,8) . '_\1 \2,\' ' . $scaff_path ;
	`$cmd`;
}


sub _get_assembly_info {
	my ($self,$file)=@_;
  open(LOG,$file) || $logger->logdie('Spades log file ' . $file . ' not found.');
  my $hash;
  while(my $line = <LOG>){
    chomp $line;
    if($line =~ /err code: -6/){
      $logger->warn('Spades assembly failed. ');
      $logger->warn($line);
      return undef;
    }
    if($line =~ / * Assembled scaffolds are in (.*)/){
      $hash->{scaffolds} = $1;
    }
  }
  close LOG;
	if(!defined($hash->{scaffolds})){
		return undef;
	}
	else{
		return $hash;
	}
}


sub _parse_spades_exec_log {
	my ($self,$res)=@_;
	my @lines = split("\n",$res);
	foreach my $line (@lines){
		if($line =~ /err code: -6/){
			return 1;
		}
	}
	return 0;
}


sub _get_unknown_reads {
	my ($self)=@_;
	my @read_list;
	foreach my $scaffold_id (@{$self->{_unknown_scaff}}){
		if(defined($self->{_blat_mapping}->{$scaffold_id})){
			# $logger->info($scaffold_id . ' containing ' .scalar(@{$self->{_blat_mapping}->{$scaffold_id}}) . ' reads.');
			push(@read_list,@{$self->{_blat_mapping}->{$scaffold_id}});
		}
		else{
			# $logger->warn('no reads for scaffold ' . $scaffold_id);
		}
	}
	my ($pairs,$singletons) = _sort_pairs_and_singletons($self,\@read_list);
	return ($pairs,$singletons);
}


sub _get_rank_reads {
	my ($self,$rank)=@_;
	my @read_list;
	foreach my $scaffold_id (@{$self->{_scaffold_by_taxo_rank}->{$rank}}){
		$logger->debug($scaffold_id);
		if(defined($self->{_blat_mapping}->{$scaffold_id})){
			$logger->debug($scaffold_id . ' containing ' .scalar(@{$self->{_blat_mapping}->{$scaffold_id}}) . ' reads.');
			push(@read_list,@{$self->{_blat_mapping}->{$scaffold_id}});
		}
		else{
			$logger->debug('no reads for scaffold ' . $scaffold_id);
		}
	}
	my ($pairs,$singletons) = _sort_pairs_and_singletons($self,\@read_list);
	return ($pairs,$singletons);
}


sub _print_pairs_and_singletons {
	my ($self,$rank,$pairs,$singletons)=@_;
	my $rank_path = cwd() . '/' . $out_path . '/' . $rank;
	if(! -e $rank_path){
		`mkdir -p $rank_path`;
	}

	open(PAIR,">$rank_path/$rank\_pairs.fa");
	open(SING,">$rank_path/$rank\_singletons.fa");

	foreach my $id (@{$pairs}){
		print PAIR $self->{_reads_index}->retrieveFastaBlock($id);
	}
	if(defined($singletons)){
		foreach my $id (@{$singletons}){
			print SING $self->{_reads_index}->retrieveFastaBlock($id);
		}
	}
	close PAIR;
	close SING;
}


sub _build_assembly_cmd {
	my ($self,$rank)=@_;
	my $rank_path = cwd() . '/' . $out_path . '/' . $rank;
	my $cmd = 'metaspades.py';
	if(! -z "$rank_path/$rank\_pairs.fa"){
		$cmd .= ' --hqmp1-12 ' . $rank_path . '/' . $rank . '_pairs.fa';
	}
	if(! -z "$rank_path/$rank\_singletons.fa") {
		$cmd .= ' --s1 ' . $rank_path . '/' . $rank . '_singletons.fa';
	}
  $cmd .= ' -o ' . $rank_path . '/' . 'spades' ;
	$cmd .= ' -t 10';
	$cmd .= ' --only-assembler --careful --cov-cutoff off';
	$self->{_cmd}->{$rank} = $cmd;
}


sub _print_reads {
	my ($self,$rank)=@_;
	# my ($pairs_unknown,$singletons_unknow) = _get_unknown_reads($self);
	my ($pairs_rank,$singletons_rank) = _get_rank_reads($self,$rank);
	# _print_pairs_and_singletons($self,$rank,$pairs_unknown,$singletons_unknow,1);
	_print_pairs_and_singletons($self,$rank,$pairs_rank,$singletons_rank,0);
	_build_assembly_cmd($self,$rank);
}


sub _sort_pairs_and_singletons {
	my ($self,$read_list)=@_;
	my $already_seen;
	my @pairs;
	my @singletons;
	foreach my $id_1 (@{$read_list}){
		my $base_1 = $1 if $id_1 =~ /(.*)\/[12]/;
		my $found=0;
		if(defined($already_seen->{$id_1})){
			next;
		}
		foreach my $id_2 (@{$read_list}){
			if($id_1 eq $id_2){
				next;
			}
			my $base_2 = $1 if $id_2 =~ /(.*)\/[12]/;
			if($base_1 eq $base_2){
				$found =1;
				if($id_1 =~ /.*\/1$/){
					push(@pairs,$id_1);
					push(@pairs,$id_2);
				}
				else{
					push(@pairs,$id_2);
					push(@pairs,$id_1);
				}
				$already_seen->{$id_2}++;
			}
		}
		if($found ==0){
			my $id_2;
			my $base_1 = $1 if $id_1 =~ /(.*)\/[12]/;
			if($id_1 =~ /.*\/1$/){
				$id_2 = $base_1 . '/2';
			}
			else{
				$id_2 = $base_1 . '/1';
			}
			if(defined($self->{_reads_index}->{index}->{$id_2})){
				if($id_1 =~ /.*\/1$/){
					push(@pairs,$id_1);
					push(@pairs,$id_2);
				}
				else{
					push(@pairs,$id_2);
					push(@pairs,$id_1);
				}
			}
			else{
				push(@singletons,$id_1);
			}
		}
		$already_seen->{$id_1}++;
	}
	$logger->info(scalar(@pairs) . ' pairs and ' . scalar(@singletons) . ' singletons.');
	return (\@pairs,\@singletons);
}


sub _sort_by_rank {
	my ($self,$desired_taxo_rank,$seek_up_or_down)=@_;
  my $rank=0;
	foreach my $hit (@{$self->{_blast_annotation}}){
    if(! defined ($hit->{taxonomy})){
      push(@{$self->{_scaffold_by_taxo_rank}->{'unknown'}},$hit->{query_id});
      next;
    }
    if($hit->{taxonomy} !~ /Virus/){
      next;
    }
		my @taxo = split(';',$hit->{taxonomy});
    $rank = $desired_taxo_rank;
		if(!defined($taxo[$rank])){
			$logger->warn('No taxo rank defined for ' . $hit->{taxonomy});
      $rank = scalar(@taxo)-1;
		}
		if($taxo[$rank] eq ''){
      if($seek_up_or_down==1){
        for(my $i = $rank+1; $i <= $#taxo ; $i++){
  				if(defined($taxo[$i])){
  					if($taxo[$i] ne ''){
  						my $rank_wo_spaces = $taxo[$i] =~ s/ /_/gr;
  						push(@{$self->{_scaffold_by_taxo_rank}->{$rank_wo_spaces}},$hit->{query_id});
              $self->{rank_to_full_taxo}->{$rank_wo_spaces} = $hit->{taxonomy};
  						last;
  					}
  				}
  				else{
  					push(@{$self->{_scaffold_by_taxo_rank}->{'no_rank'}},$hit->{query_id});
  				}
  			}
      }
      else{
        for(my $i = $rank+1; $i <= $#taxo ; $i--){
  				if(defined($taxo[$i])){
  					if($taxo[$i] ne ''){
  						my $rank_wo_spaces = $taxo[$i] =~ s/ /_/gr;
  						push(@{$self->{_scaffold_by_taxo_rank}->{$rank_wo_spaces}},$hit->{query_id});
              $self->{rank_to_full_taxo}->{$rank_wo_spaces} = $hit->{taxonomy};
  						last;
  					}
  				}
  				else{
  					push(@{$self->{_scaffold_by_taxo_rank}->{'no_rank'}},$hit->{query_id});
  				}
  			}
      }
		}
		else{
			my $rank_wo_spaces = $taxo[$rank] =~ s/ /_/gr;
			push(@{$self->{_scaffold_by_taxo_rank}->{$rank_wo_spaces}},$hit->{query_id});
      $self->{rank_to_full_taxo}->{$rank_wo_spaces} = $hit->{taxonomy};
		}
	}
  open(STATS,">abt_stats.txt");
  foreach my $k (keys(%{$self->{_scaffold_by_taxo_rank}})){
    print STATS scalar(@{$self->{_scaffold_by_taxo_rank}->{$k}}) . "\t" . $k;
    if(defined($self->{rank_to_full_taxo}->{$k})){
      print STATS "\t" . $self->{rank_to_full_taxo}->{$k};
    }
    print STATS "\n";
  }
  close STATS;
}


sub _set_options {
	my ($self)=@_;

	$self->{_selected_taxo_rank} = $selected_taxo_rank;

	if(-e $blast_ecsv_file){
		$self->{_blast_ecsv_file} = abs_path($blast_ecsv_file);
	}
	else{
		$logger->error('You must provide a ecsv blast annotation file.');
		exit;
	}
	if($bam_file eq ''){
		$logger->error('You must provide one mapping file, either psl or bam.');
		exit;
	}
	else{
		if($bam_file ne ''){
			if(-e $bam_file){
				$self->{_bam_file} = abs_path($bam_file);
			}
			else{
				$logger->error('Bam file provided not found.');
				exit;
			}
			if($pair_R1 eq '' || $pair_R2 eq ''){
				$logger->error('With Bam file, you must provide fastq files.');
				exit;
			}
			else{
				$self->{_pair_R1} = abs_path($pair_R1);
				$self->{_pair_R2} = abs_path($pair_R2);
			}
		}
	}
  if($singletons ne ''){
    $self->{_singletons} = abs_path($singletons);
  }

  if($seek_up_or_down !~ /^up$|^down$/){
    $logger->error('-seek option only accepts up|down words ');
    exit;
  }
  else{
    if($seek_up_or_down eq 'up'){
      $self->{_seek_up_or_down}=0;
    }
    else{
      $self->{_seek_up_or_down}=1;
    }
  }

  if($mega_merge){
    $self->{_mega_merging}=1;
  }
  else{
    $self->{_mega_merging}=0;
  }

	if($sample_id eq ''){
		my @a;
		if(defined($self->{_blat_psl_file})){
			@a = split('_',$self->{_blat_psl_file});
		}
		if(defined($self->{_bam_file})){
			@a = split('_',$self->{_bam_file});
		}
		$self->{_sample_id} = $a[0];
	}
	else{
		$self->{_sample_id} = $sample_id;
	}

	if($output ne ''){
		$self->{_output} = $output;
	}
	else{
		$self->{_output} = $self->{_sample_id} . '_scaffold.abt.fna'
	}
  $self->{_n_cpu} = $n_cpu;
  $self->{_read_norm} = $read_norm;
}


sub help {
my $prog = basename($0) ;
print STDERR <<EOF ;
### $prog $VERSION ###
#
# AUTHOR:     Sebastien THEIL
# VERSION:    $VERSION
# LAST MODIF: $lastmodif
# PURPOSE:    Based on annotation provided, this script will split input reads by viral family and launch assembly.
#

USAGE: perl $prog -e ecsv_file -b bam_file -1 reads.r1.fq -2 reads.r2.fq[OPTIONS]

          -e|ecsv      [EXTENDED BLAST FILE] At least two fields: query_id & taxonomy.
          -b|bam       [BAM FILE] Alignment of reads on scaffolds.
          -1           [FASTQ FILE] Pair R1.
          -2           [FASTQ FILE] Pair R2.
          ### OPTIONS ###

          -singletons  [FASTQ FILE] Single reads.
          -rn          Read normalization method.
          -seek        [up|down] If selected rank has no clade, seek for valid one up (superkingdom) or down (specie).
          -mm|megamerge          Use MegaMerge tool to merge all assemblies.
          -n_cpu       [INT] Number of CPU to use in Spades.
          -v|verbosity [INTEGER]
          -help|h      Print this help and exit.

EOF
exit(1);
}
