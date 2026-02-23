#!/usr/bin/env perl
#
#====================================================================================
#
#   Interface to run NAMD within an ONIOM QM:QM' schemme 
#
#   Developers: Federico J. Hernandez and Rachel Crespo Otero
#
#====================================================================================
#

use strict;
# use warnings;
use File::Copy;
use File::Path;
use lib join('/',$ENV{"NX"},"lib");
use colib_perl;

# We are using "strict." Thus, all variables mus be declared.
my ($mld,$mdle,$BASEDIR,$DEBUG,$JAD,$fpar,$epotf);
my ($lvprt,$nis,$nfs,$mxns,$ns,$ia);
my (@symb,@zn,@x,@y,@z,@Mass);
my (@epot,@oos);
my ($prog,%progconf);
my ($retval,$line,@line);
my ($fro_inp,$fro_soft,$fro_type,$fro_energies,$fro_oos);

# Define basic variables
define_variables();

# Read NX parameters
load_cpar();

# Read current geometry
read_geom();

# Read parameters for the TPP
read_fpar();

# Prepare TPP input
prepare_input();

# Execute TPP
run_program();

# Save files
save_files();

# Read energies from TPP output
read_energy();

# Read oscillator strengths from TPP output
read_oos();

#
#====================================================================================
#
# 			START SUBROUTINES
#
#====================================================================================
#

sub define_variables{
#
#====================================================================================
#
# Define basic variables to be used everywhere.
#
#------------------------------------------------------------------------------------
#

  $mld    = $ENV{"NX"};
  $mdle   = "run_fromage_initcond.pl:";

  $BASEDIR=`pwd`;
  chomp ($BASEDIR);

  $DEBUG  = "DEBUG";
  $JAD    = "JOB_AD";

  $prog     = getkeyword("initqp_input","prog","");
  %progconf = prog_config($prog);
  $fpar = $progconf{parfile};
  $fro_inp = "fromage.in";
  $fro_energies = "fro_energies.dat";
  $fro_oos = "fro_oos.dat";

  if (-e "zero") {
    $epotf="epot0";
  } else {
    $epotf="epot";
  }  

}

sub load_cpar{
#
#====================================================================================
#
#------------------------------------------------------------------------------------
#

  $nis     = getkeyword("initqp_input","nis" ,"");
  $nfs     = getkeyword("initqp_input","nfs" ,"");
  $lvprt   = getkeyword("initqp_input","lvprt","");

  # Check maximum between $nis and $nfs
  $mxns = max_ns($nis,$nfs);

}

sub read_geom{
#
#====================================================================================
#
#------------------------------------------------------------------------------------
#

  if ($lvprt>=3) {print_STDOUT("$mdle Reading current geometry.\n");}

  open(GM,"geom") or die "$mdle Cannot open geom";
  $ia = 0;
  while(<GM>){
    chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
    ($symb[$ia],$zn[$ia],$x[$ia],$y[$ia],$z[$ia],$Mass[$ia])=split(/\s+/,$_);
    $ia++;
  }
  close(GM);

}

sub read_fpar{}

sub prepare_input{
#
#====================================================================================
#
#------------------------------------------------------------------------------------
#

  my $au2ang = units("au2ang");

  if ($lvprt>=3) {print_STDOUT("$mdle Preparing input for third-party program.\n");}

  # Prepare input for single-point run (energy + oscillator strengths)
  
  # Transform the coordinates from NX format to XYZ format
  $retval = callprog("","$mld/nx2xyz",$mdle);
  system("mv geom.xyz mol.init.xyz");
  if ($retval != 0){
    die "$mdle is dying now\n";

}

  my $fro_soft = getfro_soft($fro_inp);
  # Check if $fro_soft is empty (high_level not found)  
  if ($fro_soft eq "") {
    print_STDOUT("fromage selected software is Gaussian \n");
  }else {
    print_STDOUT("fromage selected software is $fro_soft \n");
  }

}

sub getfro_soft{
#
#====================================================================================
# This subroutine recognizes the electronic structure software used by fromage to
# model the QM region. It saves the software name to later call the proper
# run_cioverlap code. For now it is only implemented for Gaussian and Turbomole.
#
# -----------------------------------------------------------------------------------
  # Initialize $fro_soft as empty string 
  my $fro_soft = "";

  # Open the file for reading
  open(IN,$fro_inp) or die "Cannot open fromage.in";

  # Read the file line by line
  my(@line);
  while (<IN>) {
    chomp $line;  # Remove trailing newline character
    @line = split(/\s+/,$_); 
    if ($line[0] = "high_level") {
      # Extract the value in the second column and store it in $fro_soft
      $fro_soft = $line[1];
      last;  # Exit the loop since we found the desired line
    }
  }

  close (IN);
  return $fro_soft;
}

sub get_fro_type{
#
#====================================================================================
#
# This subroutine recognizes the type of ONIOM calculation. If the nuclei in the QM'
# region will remain clamped or will allow to evolve. Check fromage documentaiton
# to see how to set a flexible QM' region.
#
# -----------------------------------------------------------------------------------
  my $fro_type = "";
  my ($nat_flex);

  # Open the file for reading
  open(IN,$fro_inp) or die "Cannot open fromage.in";

  # Read the file line by line
  $nat_flex=0;
  my(@line);
  while (<IN>) {
    chomp $line;  
    @line = split(/\s+/,$_); 
    if ($line[0] = "ll_flex_natoms") {
      $nat_flex = $line[1];
      last;
    }
  }
  close (IN);

  if ($nat_flex == 0) {
    $fro_type = "fro_run.py";
  }else {
    $fro_type = "fro_run_flex.py";
  }

  return $fro_type;

}

sub run_program{
#
#====================================================================================
#
#------------------------------------------------------------------------------------
#
	
  my ($retval);

  if ($lvprt>=3) {print_STDOUT("$mdle Executing fromage program.\n");}

#  my $tpp = $ENV{"My_TPP"};

  my $fro_type = get_fro_type($fro_inp);

  $retval = callprog("",$fro_type,$mdle);
  print_STDOUT("$mdle fromage selected $fro_type.\n")
#  if ($retval != 0){print_STDOUT("Error in fromage program execution!\n");}

  # check_for_problems();

}

sub read_energy{
#
#====================================================================================
#
#------------------------------------------------------------------------------------
#

  if ($lvprt>=3) {print_STDOUT("$mdle Reading potential energies.\n");}

  # This is just a pseudocode suggestion. It must be adapted to the actual case.
  open(IN,$fro_energies) or die "Cannot find fro_energies.dat";
  $ns = 0;
  my(@line);
  while(<IN>){
    chomp($_);
    $_ =~ s/^\s+//;
    @line = split(/\s+/,$_);
    $epot[$ns] = $line[0];
    $ns++;
  }

  close(IN);

  write_energy();

}

sub write_energy{
#
#====================================================================================
#
#------------------------------------------------------------------------------------
#

  my ($ifs);

  if ($lvprt>=3) {print_STDOUT("$mdle Writing potential energies.\n");}

  for ($ns = 1; $ns <= $nfs-$nis; $ns++){
    $ifs = $nis + $ns;
    open(OUT,">$epotf.$ifs") or die "Cannot write to $epotf.$ifs";
    print OUT "$epot[0]\n";
    print OUT "$epot[$ns]\n";
    close(OUT);
  }

}

sub read_oos{
#====================================================================================
#
#------------------------------------------------------------------------------------
#

  if ($lvprt>=3) {print_STDOUT("$mdle Reading oescillator strengths.\n");}

  # This is just a pseudocode suggestion. It must be adapted to the actual case.
  open(IN, $fro_oos) or die "Cannot find fro_oos.dat output";
  $ns = 0;
  my(@line);
  while(<IN>){
    chomp($_);
    $_ =~ s/^\s+//;
    @line = split(/\s+/,$_);
    $oos[$ns] = $line[0];
    $ns++;
  }

  close(IN);

  write_oos();

}

sub write_oos{
#
#====================================================================================
#
#------------------------------------------------------------------------------------
#
  my ($ifs);

  for ($ns = 0; $ns <= $nfs-$nis-1; $ns++){
    $ifs = $nis + $ns + 1;
    open(OUT,">oos.$ifs") or die "Cannot write to oos.$ifs";
    print OUT "$oos[$ns]\n";
    close(OUT);
  }

}
