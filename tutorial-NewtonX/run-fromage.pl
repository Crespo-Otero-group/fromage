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
my ($mld,$mdle,$BASEDIR,$DEBG,$JAD,$JND,$ctd,$fpar);
my ($mem,$nat,$istep,$nstat,$nstatdyn,$ndamp,$kt,$dt,$t,$tmax,$nintc);
my ($nxrestart,$thres,$killstat,$timekill,$prog,$lvprt,$etot_jump,$etot_drift);
my (@symb,@zn,@x,@y,@z,@Mass,$ia);
my ($typeofinput,$type,$jobtype);
my (@epot,$ns);
my ($vdoth,$run_flag);
my ($ncoup,$nc);
my (@gx,@gy,@gz,@hx,@hy,@hz);
my ($kross_d,$cascade_d,$current_d,$never_state_d,$include_pair_d);
my ($e_ci_d,$ci_cons_d,$cio_options_d,$cisc_options_d,$idalton_d);
my ($cprog_d,$coptda_d,$ncore_d,$ndisc_d,$blasthread_d);
my ($blasthread,$CIO_OPTIONS,$CISC_OPTIONS,$CPROG);
my ($retval,$file,$epot);
my (%progconf);
my ($fro_inp,$fro_energies,$fro_gradients,$fro_nacs,$fro_socs,$fro_soft,$fro_type,$qm_natoms);
my ($line,$values,$keyword1_d,$keyword2_d,$keyword1,$keyword2,$IN,$js,@line);
my ($mldcio,$from,$RS,$eps);

# Define basic variables
define_variables();

# Write debug messages
write_debug1();

# Read current geometry
read_geom();

# Read TPP parameters
read_fpar();  

# Prepare TPP input
prepare_input();

# Run TPP
run_program();

# Read and write energies
read_energy();

# Read and write gradients
read_gradients();

# Nonadiabatic couplings
treat_nacme();

#
#====================================================================================
#
# 				START SUBROUTINES
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
  $mdle   = "run-fromage.pl:";  # %% CHANGE HERE %%

  $BASEDIR=`pwd`;
  chomp ($BASEDIR);

  $RS    = "RESULTS";
  $DEBG  = "DEBUG";
  $JAD    = "JOB_AD";
  $JND    = "JOB_NAD";
  $ctd    = "control.d";
  $fro_inp = "fromage.in";
  $fro_energies = "fro_energies.dat";
  $fro_gradients = "fro_gradients.dat";
  $fro_nacs = "fro_nacs.dat";
  $fro_socs = "fro_socs.dat";

  $eps = 1E-9;

  # Load current dynamics status (read control.d)
  ($nat,$istep,$nstat,$nstatdyn,$ndamp,$kt,$dt,$t,$tmax,$nintc,
   $mem,$nxrestart,$thres,$killstat,$timekill,$prog,$lvprt,
   $etot_jump,$etot_drift)=load_status($ctd,$mdle); 

  %progconf = prog_config($prog);
  $fpar = $progconf{parfile};  
}

sub write_debug1{
#
#====================================================================================
#
# Write debug messages.
#
#------------------------------------------------------------------------------------
#
  if ($lvprt>=3) {print_STDOUT("$mdle has taken over\n",$istep,$kt);}
  if ($lvprt>=3) {print_STDOUT("$mdle running in $BASEDIR\n",$istep,$kt);}
  if ($lvprt>=3) {print_STDOUT("$mdle beginning here\n",$istep,$kt);}
  if ($lvprt>=3) {
     print_STDOUT("$mdle Dynamics control: \n",$istep,$kt);
     print_STDOUT("$nat $istep $nstat $nstatdyn $ndamp $kt $dt $t $tmax \n",$istep,$kt);
     print_STDOUT("$nintc $mem $nxrestart $thres $killstat $timekill $prog $lvprt\n",$istep,$kt);
     print_STDOUT("$etot_jump $etot_drift\n\n",$istep,$kt);
  }
}

sub read_geom{
#
#====================================================================================
#
#------------------------------------------------------------------------------------
#

  if ($lvprt>=3) {print_STDOUT("$mdle Reading current geometry.\n",$istep,$kt);}

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

sub fro_get_nacs{
#
#====================================================================================
#
# This subroutine reads and organizes the nacs when they are explicitly computed 
#
#------------------------------------------------------------------------------------

  my($nactype,$coupling,@nacs,$j,@nac_pair,$first,$second,$out_file);
  $j=0;
  open(IN, "transmomin") || die "Cannot open transmomin\n";
  while(<IN>) {
    if (!/CI/) {
      $j++;
      chomp($_);
      @nac_pair = split(/\s+/,$_); # Split line at whitespace
      push(@nacs,"$nac_pair[3] $nac_pair[1] \n");
    }
  }
  close(IN);

  # Open the file in append mode
  open(my $out_file, '>>', $fro_inp) or die "Could not open file '$fro_inp'";
  
  # Print the coupling value followed by pairs from @nac_pair
  print $out_file "$coupling ";
  foreach my $pair (@nac_pair) {
    my ($first, $second) = @$pair;  # Dereference the array reference
    print $out_file "$first $second ";
  }
  close($out_file);
}
  
sub getfro_qmnatoms{
#
#====================================================================================
#
# Get the number of atoms in the QM region to set the CIOVERLAP calculation
#
# -----------------------------------------------------------------------------------

 my ($fro_inp) = @_;
 my $qm_natoms = "";

 open(my $IN, '<', $fro_inp) or die "Cannot open $fro_inp: $!";
 
# my(@line);
  while (<$IN>) {
    chomp;  # Remove trailing newline character
    my @line = split(/\s+/);
    if ($line[0] eq "hl_natoms") {
      $qm_natoms = $line[1];
      last;
    }
  }

  close ($IN);
  return $qm_natoms;
}

sub getfro_soft{
#
#====================================================================================
#
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
    if ($line[0] eq "high_level") {
      # Extract the value in the second column and store it in $fro_soft
      $fro_soft = $line[1];
      last;  # Exit the loop since we found the desired line
    }
  }

  # Close the file
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
    if ($line[0] eq "ll_flex_natoms") {
      $nat_flex = $line[1];
      last; 
    }
  }
  close (IN);

  if ($nat_flex == 0) {
    $fro_type = "fro_run.py";
  }else {
    $fro_type = "fro_run.py"; # fro_run_flex.py
  }

  return $fro_type;
  
}  

sub prepare_input{
#====================================================================================
#
#------------------------------------------------------------------------------------
#

  my $au2ang = units("au2ang");

  if ($lvprt>=3) {print_STDOUT("$mdle Preparing input for fromage program.\n",$istep,$kt);}

  find_typeofinput();

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

  if ($typeofinput eq "jad"){
    # Prepare input for adiabatic run (energy + gradients like in JOB_AD/)
    if ($lvprt>=3) {print_STDOUT("$mdle Adiabatic dynamics with fromage.\n",$istep,$kt);}
  }elsif($typeofinput eq "jnd"){
    # Prepare input for nonadiabatic run (energy + gradients + NACME like in JOB_NAD/)
    # If TPP does not provide NACME, you do not need to enter anything here.
    $vdoth = getkeyword("sh.inp","vdoth",2);
    if ($vdoth == 0){
      print_STDOUT("fromage will provide NACs\n",$istep,$kt);
#      fro_get_nacs();
    }
  }
}

sub find_typeofinput{
#
#====================================================================================
#
# Check type of dynamics.
#
#------------------------------------------------------------------------------------
#
 
  if ($lvprt>=3) {print_STDOUT("$mdle Checking type of dynamics.\n",$istep,$kt);}
 
  # Read $type from type_of_dyn.out 
  my $typeout="type_of_dyn.out";
  $type = 2;
  if (-s $typeout){
    open(INP,$typeout) or die "$mdle Cannot open $typeout";
    $_=<INP>;
    chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
    my $type = $_;
    close(INP);
  }

  # Define type of input
  if ((-e $JAD) and (!-e $JND)){
    $typeofinput = "jad";
  }elsif((!-e $JAD) and (-e $JND)){
    $typeofinput = "jnd";
  }elsif((-e $JAD) and (-e $JND)){
    if ($type == 1){
      $typeofinput = "jad";
    }elsif($type != 1){
      $typeofinput = "jnd";
    }
  }

}

sub run_program{
#
#====================================================================================
#
#------------------------------------------------------------------------------------
#
   
  if ($lvprt>=3) {print_STDOUT("$mdle Executing fromage program.\n",$istep,$kt);}

  # my $tpp = $ENV{"My_TPP"}; 

  my $fro_type = get_fro_type($fro_inp); 

#  if ($lvprt>=3) {print_STDOUT("$mdle LLAMARA A $fro_type.\n",$istep,$kt);}

  print_STDOUT("$mdle Calling $fro_type \n");
  callprog("",$fro_type,$mdle);
#  $retval = callprog("",$fro_type,$mdle);

  write_fro_dyn();

}

sub write_fro_dyn{
#
#====================================================================================
#
# This function writes in RESULTS/ the xyz coordinates of the complete ONIOM Embedded 
# Cluster (OEC) at every step of the dynamics
#
# -----------------------------------------------------------------------------------
  my $fro_filename = 'geom_cluster.xyz';
  my $fro_output = 'fromage_dyn.xyz';
  my $fro_file_path = join '/', "../$RS", $fro_output;
  $fro_file_path =~ s/\/\//\//g;
  my @lines;

  # Open the input file in read mode
  open(my $fh, '<', $fro_filename) or die "Cannot open file $fro_filename";
  while (my $line = <$fh>) {
    push(@lines, $line);
  }
  close($fh);
  # Replace the second line with the current time information
  $lines[1] = "Time = $t\n";
  # Write the changes to the fromage_dyn.xyz file
  open(my $fh_out, '>>', $fro_file_path) or die "Cannot open file $fro_file_path";
  print $fh_out @lines;
  close($fh_out);

}

sub read_energy{
#
#====================================================================================
#
#------------------------------------------------------------------------------------
#

  if ($lvprt>=3) {print_STDOUT("$mdle Reading Epot.\n",$istep,$kt);}
	
  open(IN,$fro_energies) or die "Cannot open fro_energies.dat";
  print_STDOUT("Energies (Eh):\n",$istep,$kt);
  $ns = 0;
  my(@line);
  while (<IN>) {
    chomp($_);
    $_ =~ s/^\s+//;
    @line = split(/\s+/,$_); # Split line at whitespace
    $epot[$ns] = $line[0];
    $ns++;
    print_STDOUT("State $ns : $line[0] \n",$istep,$kt);
  }
  close(IN);

  write_energy();
}

sub write_energy{
#
#====================================================================================
#
# This subroutine writes the potential energies to epot file.
# It supposes that @epot vector contains these energies.
# See subroutine read_energies();
# oldepot and new epot files, used for interpolation are also
# written.
#
#------------------------------------------------------------------------------------
#
  
  if ($lvprt>=3) {print_STDOUT("$mdle Writing Epot.\n",$istep,$kt);}

  $file = "epot";
  if (-s $file){
    copy($file,"oldepot") or die "Copy failed: $!";
  } 
  open(OUT,">$file") or die "Cannot write to $file";
  foreach(@epot){
    print OUT "$_\n"; 
  }
  close(OUT);
  copy($file,"newepot") or die "Copy failed: $!";
}

sub read_gradients{
#
#====================================================================================
#
#------------------------------------------------------------------------------------
#

  if ($lvprt>=3) {print_STDOUT("$mdle Reading Gradients.\n",$istep,$kt);}

  # Initialize array with zeros
  for ($ns = 0; $ns <= $nstat-1; $ns++){
    for ($ia = 0; $ia <= $nat-1; $ia++){
      $gx[$ns][$ia] = 0.0;
      $gy[$ns][$ia] = 0.0;
      $gz[$ns][$ia] = 0.0;
    }
  }

  open(IN,$fro_gradients) or die "Cannot find fro_gradients.dat";
  $ns = 0;
  $ia = 0;
  my(@line);
  while (<IN>) {
    chomp($_);
    $_ =~ s/^\s+//;
    @line = split(/\s+/,$_);

    $gx[$ns][$ia] = sprintf("%.15f",$line[0]);
    $gy[$ns][$ia] = sprintf("%.15f",$line[1]);
    $gz[$ns][$ia] = sprintf("%.15f",$line[2]);

    $ia++;
    if ($ia >= $nat) {
      $ia = 0;   # Reset atom index
      $ns++;    # Move to the next state 
    }
  }
  close(IN);

  $ia = 0;
  if ($lvprt>=3) {
    print_STDOUT("Gradient of current state = $nstatdyn: \n",$istep,$kt);
    for ($ia = 0; $ia <= $nat-1; $ia++){
      print_STDOUT("$gx[$nstatdyn-1][$ia]   $gy[$nstatdyn-1][$ia]   $gz[$nstatdyn-1][$ia] \n",$istep,$kt);
    }
  }

  write_grad();

}

sub write_grad{
#
#====================================================================================
#
# This subroutine writes the energy gradients to grad and grad.all files.
#
#------------------------------------------------------------------------------------
#

  if ($lvprt>=3) {print_STDOUT("$mdle Writing Gradients.\n",$istep,$kt);}
 
  open(OUT1,">grad.all") or die "Cannot write to grad.all";
  open(OUT2,">grad")     or die "Cannot write to grad";

  for ($ns = 0; $ns <= $nstat-1; $ns++){
    for ($ia = 0; $ia <= $nat-1; $ia++){
      printf OUT1 "%20.14f %20.14f %20.14f\n",$gx[$ns][$ia],$gy[$ns][$ia],$gz[$ns][$ia];
      if ($ns == $nstatdyn-1){
        printf OUT2 "%20.14f %20.14f %20.14f\n",$gx[$ns][$ia],$gy[$ns][$ia],$gz[$ns][$ia];
      }
    }
  }

  close(OUT1);
  close(OUT2);
}

sub treat_nacme{
#
#====================================================================================
#
# This subroutine processes the NACME. 
# If it was computed by TPP, it reads and writes them.
# If not, it may call cioverlap to get the couplings.  
#
#------------------------------------------------------------------------------------
#
 
  if ($lvprt>=3) {print_STDOUT("$mdle Treating couplings.\n",$istep,$kt);}

  $run_flag=read_single_value("which_run_is_that","");

  # $run_flag = "second run" is a new program call after hopping.
 
  if ($run_flag ne "second run"){

    # ncoup is half of the number of coupling vectors
    $ncoup = $nstat*($nstat-1)/2;

    # Intialize arrays
    for ($nc = 0; $nc <= $ncoup-1; $nc++){
      for ($ia = 0; $ia <= $nat-1; $ia++){
         $hx[$nc][$ia] = 0.0;
         $hy[$nc][$ia] = 0.0;
         $hz[$nc][$ia] = 0.0;
      }
    }

    $vdoth = getkeyword("sh.inp","vdoth",2);
    if ($vdoth == 0){
       read_nacme();
       write_nacme("nad_vectors");
    }elsif(($vdoth == 1) or ($vdoth < 0)){
       run_cioverlap();
    }  # Even if you don't define any coupling here, you can still use TD-BA
    
  }
}

sub read_nacme{
#
# %% CHANGE HERE %%
#====================================================================================
#
#------------------------------------------------------------------------------------
#

  if ($lvprt>=3) {print_STDOUT("$mdle Reading coupling vectors.\n",$istep,$kt);}

  my ($is,$js,$nc,$ia);

  open(IN,$fro_nacs) or die "Cannot find fro_nacs.dat";
  $ncoup = $nstat*($nstat-1)/2;
  $nc = 0;
  $ia = 0;
  my(@line);
  while (<IN>) {
    chomp($_);
    $_ =~ s/^\s+//;
    @line = split(/\s+/,$_);

    $hx[$nc][$ia] = sprintf("%.15f",$line[0]);
    $hy[$nc][$ia] = sprintf("%.15f",$line[1]);
    $hz[$nc][$ia] = sprintf("%.15f",$line[2]);

    $ia++;
    if ($ia >= $nat) { 
      $ia = 0;   # Reset atom index
      $nc++;    # Move to the next pair of states
    }
  }
  close(IN);
}

sub write_nacme{
#
#====================================================================================
#
# This subroutine writes @h to a file. The definition and order  of @h are explained
# in subroutine @read_nacme. The output file is passed as an argument of the
# subroutine.
#
#------------------------------------------------------------------------------------
#
 
  if ($lvprt>=3) {print_STDOUT("$mdle Writing coupling vectors.\n",$istep,$kt);}

  ($file) = @_;  
  open(OUT,">$file") or die "Cannot write to $file";

  for ($nc = 0; $nc <= $ncoup-1; $nc++){
    for ($ia = 0; $ia <= $nat-1; $ia++){
      print OUT "$hx[$nc][$ia]  $hy[$nc][$ia]  $hz[$nc][$ia]\n";
    }
  } 

  close(OUT);  

}

sub run_cioverlap{
#
#====================================================================================
#
#------------------------------------------------------------------------------------
#
 
  if ($lvprt>=3) {print_STDOUT("$mdle Computing state-overlap matrix.\n",$istep,$kt);}
 
  # The main variable determining the behavior of the cioverlap 
  # programs are in jiri.inp file.
  #

  ($kross_d,$cascade_d,$current_d,$never_state_d,$include_pair_d,$e_ci_d,$ci_cons_d,
   $cio_options_d,$cisc_options_d,$idalton_d,$cprog_d,$coptda_d,
   $ncore_d,$ndisc_d,$blasthread_d)=load_defaults("jiri.inp",$prog,$vdoth);

  my $fro_soft = getfro_soft($fro_inp);
  my $qm_natoms = getfro_qmnatoms($fro_inp);

  call_cioverlap($fro_soft, $qm_natoms);
 
  # Read cioverlap results and write them to nad_vectors or use them for local diabtization
  #
  if ($vdoth == 1){
    if ($lvprt>=3) {print_STDOUT("Overlap matrix will be used to compute time-derivative couplings.\n",$istep,$kt);}
    $retval = callprogsp($mld, "read_cioverlap.pl", $mdle);
    if ($retval != 0){die "$mdle is dying now (read_cioverlap.pl)\n";}
  }elsif($vdoth < 0){
    if ($lvprt>=3) {print_STDOUT("Overlap matrix will be used to do local diabatization.\n",$istep,$kt);}
    # In this case, nothing needs to be done. sh program will automatically read run_cioverlap.log.
  }

}

sub change_control_d{
#
#====================================================================================================
# This function modifies temporarily the contrl.d file by changing the number of active atoms (those
# that evolve during the dynamics) by the number of atoms only in the QM region. This is to perform
# the time-dependent overlap
#
# ----------------------------------------------------------------------------------------------------

  my ($qm_natoms) = @_;
  
  # Modify the first entry in the original file
  open(my $fh, '+<', $ctd) or die "Failed to open $ctd for modification: $!";
  # Read the first line and replace the first number
  my $line = <$fh>;
  if ($line =~ /^(\s*\d+),/) {
     $line =~ s/^(\s*\d+),/$qm_natoms,/;
  }

  seek($fh, 0, 0) or die "Failed to rewind $ctd: $!";
  print $fh $line;
  truncate($fh, tell($fh)) or die "Failed to truncate $ctd: $!";  # Truncate at the current position

  close($fh);

}

sub truncate_geom {
#
#====================================================================================================
# This function truncates the file geom to the atoms only in the QM region. This is only use to compute
# the time-dependent couplings via the overlaps. Then the file geom is restored to its original form
#
# ----------------------------------------------------------------------------------------------------
#

  my ($qm_natoms) = @_;
  open(my $GM, '<', 'geom') or die "Cannot open geom: $!";
  my @lines = <$GM>; 
  close($GM);

  if ($qm_natoms > @lines) {
    die "Requested number of atoms (\$qm_natoms = $qm_natoms) exceeds the number of lines in 'geom'.";
  }
  
  # Keep only the first $qm_natoms lines
  my @selected_lines = @lines[0 .. $qm_natoms - 1];

  # Write the truncated content back to the 'geom' file
  open(my $OUT, '>', 'geom') or die "Cannot open geom for writing: $!";
  print $OUT @selected_lines;
  close($OUT);
  
}

sub call_cioverlap {
#
#====================================================================================================
# This function recognises the electronic structure software used by fromage to describe
# the QM region, organises all the necesary information and calls the corresponding NX
# run_cioverlap_$software.pl function
# 
#----------------------------------------------------------------------------------------------------
#
  my ($fro_soft, $qm_natoms) = @_;  # Accept two arguments

  system("cp control.d control.d_fro_tmp");
  system("cp geom geom_fro_tmp");
  change_control_d($qm_natoms);
  truncate_geom($qm_natoms);
  
  if (-s "tmp.old/WORK/phases" ){
    system("cp -f tmp.old/WORK/phases WORK/phases.old");
  }

  if ($fro_soft eq 'gaussian'){
    if ($istep == 0){
      inputgau(); # calling intputgau from run-gau.pl
                  # to extract from the gaussian input file
                  # information required from run_cioverlap_gau.pl
    }
   
    read_gaussian_info();

    if ($lvprt>=3) {print_STDOUT("$mdle Calling run_cioverlap_gau.pl.\n",$istep,$kt);}

    $retval = callprog($mld,"run_cioverlap_gau.pl > run_cioverlap.out",$mdle);
    if ($retval != 0){
      die "$mdle Error in run_cioverlap_gau.pl\n";
    }
    if ($lvprt>=2){
      system("cp -f run_cioverlap.log $BASEDIR/../$DEBG/run_cioverlap.log.$t");
      system("cp -f run_cioverlap.out $BASEDIR/../$DEBG/run_cioverlap.out.$t");
    }
    system("mv control.d_fro_tmp control.d"); # Restores the original form of control.d 
    system("mv geom_fro_tmp geom"); # Restores the original form of geom
    print_STDOUT("\n $mdle //// NONAD \n",$istep,$kt);
    read_write_nonad();

  }elsif (($fro_soft eq 'turbomole') or ($fro_soft eq 'turbomole_tddft')) {
    system("mv ricc2.out grad.out");
    $retval = callprog($mld,"run_cioverlap_turbo.pl > run_cioverlap.out",$mdle);
    if ($retval != 0){
      die "$mdle Error in run_cioverlap_turbo.pl\n";
    }
    if ($lvprt>=2){
      system("cp -f run_cioverlap.log $BASEDIR/../$DEBG/run_cioverlap.log.$t");
    }
    system("mv control.d_fro_tmp control.d"); # Restores the original form of control.d
    system("mv geom_fro_tmp geom"); # Restores the original form of geom
    print_STDOUT("\n $mdle //// NONAD \n",$istep,$kt);

    read_turbo_info();

  }elsif ($fro_soft eq 'dftb')  {
    $retval = callprog($mld,"run_cioverlap_dftb+.pl > run_cioverlap.out",$mdle);
    if ($retval != 0){
      die "$mdle Error in run_cioverlap_dftb+.pl\n";
    }
    if ($lvprt>=2){
      system("cp -f run_cioverlap.log $BASEDIR/../$DEBG/run_cioverlap.log.$t");
    }
    system("mv control.d_fro_tmp control.d"); # Restores the original form of control.d
    system("mv geom_fro_tmp geom"); # Restores the original form of geom
    print_STDOUT("\n $mdle //// NONAD \n",$istep,$kt);
    read_write_nonad();

  }elsif ($fro_soft eq 'molcas') {
     die "run_cioverlap_molcas.pl is not implemented as standalone function yet\n";
  }elsif ($fro_soft eq 'orca') {
     die "run_cioverlap_orca.pl is not implemented as standalone function yet\n";
  }else {
    die "The method $fro_soft has not been implemented yet for NX+fromage dynamics \n";
  }
}

#=====================================================================================
#============================== SOME NX SUBROUTINES ==================================
#
#------------------------ Routines for fromage + Gaussian ----------------------------
#
sub inputgau{
#
  my ($gaucom);
  $gaucom="gaussian.com";
  my($blank,$bh,$ib,$ib2,$gausskey4,$gi,$gausskey3,$dch);
  my(@g,@g2,@g5,@g6,@gw,@gausskey2,@gausskey1,@line,@ing);
  my (@gm,$defM,$method,$bst);
  $bh=0;
  $ib=0;
  $ib2=0;
  $blank=0;

  open(GC,"$gaucom") or die "$mdle Cannot open $gaucom";

  open (GKEY,">>gausskeywords") or die "Cannot open gausskeywords";
  open (FINP1,">>gaussinp1") or die "Cannot open gaussinp1";
  open (FINP2,">>gaussinp2") or die "Cannot open gaussinp2";
  open (MEM,">>memo") or die "Cannot open memo";
  open(IMu,">multiplicity.dat");  # initial multiplicity provided by the user

  print FINP1 "%chk=gaussian.chk\n";
  print FINP1 "%rwf=gaussian.rwf\n";

  while(<GC>){
    chomp;
    $_ =~ s/^\s*//; # remove leading blanks
    $_ =~ s/\s*$//; # remove trailing blanks

    if ($_ eq ""){
      $blank++;
    }

    if ($blank>=1 and $blank<=2 and $bh<=1 ) {
      @gm=split(/\s+/,$_);
      print FINP2 "$_\n";
      $defM=$gm[1];   #multiplicity given by the user   
    }
    if (/\%mem/i) {    # to avoid memory problems in the double molecule calculation
      print FINP1 "$_\n";
      print MEM "$_\n";
      close (MEM);
    }
    if (/\%nproc/i) {
      print FINP1 "$_\n";
    }

    if(/\#/){
      @g=split(/\s+/,$_);
      foreach $gi (@g){
      if ($gi!~/TD/i and $gi!~/CIS/i and $gi!~/force/i and $gi!~/\#/ and $gi!~/GEN/i and $gi!~/Nosym/i){
        push(@gausskey1,$gi);     #here the keywords can be checked
      }elsif ($gi=~/TD/i or $gi=~/CIS/i){
        @g5=split(/\(/,$gi);
        @g6=split(/\)/,$g5[1]);
        @g2=split(/\,/,$g6[0]);
        $method=$g5[0];
        open(MT,">method");
        print MT $method;
        close(MT);
        foreach (@g2){
          if ($g2[$ib]!~/Nstates/i and $g2[$ib]!~/root/i){
            $gausskey3=$g2[$ib];
            push(@gausskey2,$gausskey3);
            $ib2++;
          }
          $ib++;
        }
      }elsif ($gi=~/GEN/i){
        push(@gausskey1,$gi);
        open(BS, "basis") or die "$mdle Cannot open $bst!\nsystems answer was: $!\n";
        open (BS2,">basis2") or die "can't open basis2";
        open (DB,">>dbasis") or die "can't open dbasis"; # only need it for the double molecule calc
        while (<BS>){
          chomp $_;
          $_ =~ s/^\s*//;
          $_ =~ s/\s*$//;
          if ($_ !~ /^$/){       #print only the non-empty lines
            print BS2 "$_\n";
            print DB "$_\n";
            push(@line,$_);
          }
        }
        @ing=split(/\s+/,$line[0]);

        if ($ing[0]=~ /^[+-]?\d+$/) {
          foreach (@line){
            @gw=split(/\s+/,$_);
            if ($gw[0]=~ /^[+-]?\d+$/){
              $dch=$gw[0]+$nat;
              print DB "$dch  0\n";
            }else{
            print DB "$_\n";
            }
          }
        }
        print BS2 "           \n";
        print DB "           \n";
        }
        if ($ib2>1) {
          $gausskey4=join(',',@gausskey2);
        }else{
        $gausskey4=$gausskey2[0];
        }
      }
    }
    if ($blank eq 2){
      $bh++;
    }
  }

  print GKEY "@gausskey1\n";
  print IMu "$defM\n";

  if ($gausskey4) {
    print GKEY "$gausskey4\n";
  }

  close(GC);
  close(GKEY);
  close(FINP1);
  close(FINP2);
  close(MEM);
  close(BS);
  close(BS2);
  close(DB);
  close(IMu);
}


sub read_gaussian_info{
#
#=============================================================================================
# Routine adapted by Federico J. Hernandez from read_write_energy_grad()
# found in from run-gau.pl
# Check whether the calculation converged and stop (or warn) if necessary.
# Write the main configuration (character of the state, orbital, etc), 
# if available to the NX log file.
# Create files to the interpolation routine (only for nonadiabatic cases).
#
#---------------------------------------------------------------------------------------------
#
  my (@g,$ie,$ih,$ia,$itop);
  my ($sign,$ndown,$nup,$gaulog,$kind_gau,$RS,$jobtype,@line,$energy0,$S2S0,$MSS0);
  my (@Fs,$nstatdyn_ex,$Fscurr,$energy,$au2ev,@mult,@s2,$mseek,$au2ang,$MS02);
  my ($MS0,$S0);

  $gaulog="gaussian.log";
  $au2ev = units("au2ev");
  $au2ang = units("au2ang");
  $ie=0;
  $ih=0;
  $itop=0;
  $sign=0;

  $ndown=0;
  $nup=0;

  open(GL,$gaulog) or die "$mdle Cannot open $gaulog";
  open(EN,">energies") or die "$mdle Cannot write energies.";  # energies file to be used in run-cio
  open(FLOG,">>../$DEBG/log.conv");

  # Reading Multipliticity M=2S+1
  open(IMu,"multiplicity.dat");
  $MS0=<IMu>;  #desired multiplicity
  if ($MS0==1) {
    $kind_gau=0; # HF for closed shell
  }else{
    $kind_gau=1; # UHF for open shell systems
  }

  $S0=($MS0-1)/2;  #desired spin  
  $MS02=$S0*($S0+1); #desired <S**2>

  if ($kind_gau==1){
    open(TSP,">>../$RS/time-spin.dat");  # file with S and <S**2>
    open(TSF,">>../$RS/time-F.dat");
    print TSP " $t  ";
    print TSF " $t   ";
  }

  if ($thres > 0){
     $jobtype="nonadiabatic";
  }

  if ($jobtype eq "nonadiabatic" and $istep==0){
   open(NADA,">nada");
  }

  print EN "$nstat\n";
  while(<GL>){
    if (/Convergence criterion not met/){
      print FLOG "DFT convergence: NO\n";
      exit_error("Convergence failure in the DFT calculation",$mdle);
      die;
      $ih==1
    }
    if (/SCF Done:/ and $ih==0 ){
      chomp;
      $_ =~ s/^\s*//;
      $_ =~ s/\s*$//;
      @line=split(/\s+/,$_);
      $energy0=$line[4];
      print EN "$energy0 ";
      print FLOG "DFT convergence: YES\n";
    }elsif ($kind_gau==1 and /Sz/ and /\<S\*\*2>/ and not /Initial/){   # <S**2> ground state, just if it's neeeded
      $_ =~ s/^\s*//;
      $_ =~ s/\s*$//;
      @g=split(/\s+/,$_);
      $S2S0=$g[7];
      $MSS0=2*$g[9]+1;
      $Fs[0]=($S2S0-$MS02)/(2*$MS02+2);  #Fs for the ground state
      print TSP  "$MSS0   $S2S0  "; # GS S2
      printf TSF  ("% 18.12f"), $Fs[0];
      if  ($nstatdyn_ex==0) {
        $Fscurr=$Fs[0];
      }
      $ih=0;
    }elsif (/Excited State /){
      chomp;
      $_ =~ s/^\s*//;
      $_ =~ s/\s*$//;
      @line=split(/\s+/,$_);
      $energy=$energy0+$line[4]/$au2ev;
      print_STDOUT("$_\n",$istep,$kt);
      print EN "$energy ";

      while ($itop == 0){    # writing contributions
        $_=<GL>;
        if(/\->/ or /\<-/ ) {
          print_STDOUT("$_",$istep,$kt);
        }else {
          $itop=1;
        }
      }
      $itop=0;

      $ie++;
      print FLOG "TDDFT State $ie : YES\n";
      if ($kind_gau==1) {  # New cycle
        @mult=split(/-/,$line[3]);  # Multiplicity
        @s2=split(/=/,$line[9]);     # <S**2>
        print_STDOUT("Multiplicity $mult[0]  \<S\**2>  $s2[1] \n");
        print TSP  "$mult[0]     $s2[1]     ";
        $Fs[$ie]=($s2[1]-$MS02)/(2*$MS02+2);  #Fs for the excited states
        printf TSF  ("% 18.12f"), $Fs[$ie];
        if  ($nstatdyn_ex==$ie) {
          $Fscurr=$Fs[$ie];
        }
      }
    }
    if ($jobtype eq "nonadiabatic" and $istep==0){ # reading for nonadiabatic dynamics only in the first step
      if (/NBasis/ and /NAE/ and /NBE/ and /NFC/ and /NFV/){
        $_ =~ s/^\s*//;
        $_ =~ s/\s*$//;
        @g=split(/\s+/,$_);
        print NADA "$g[1] "#$Nbf
      }elsif (/NROrb/ and /NOA/ and /NOB/ and /NVA/){
        $_ =~ s/^\s*//;
        $_ =~ s/\s*$//;
        @g=split(/\s+/,$_);
        print NADA "$g[3] $g[5] $g[7] $g[9] "; #"$NoccA $NoccB $NvirtA $NvirtB";
      }elsif (/roots to seek/){
        $_ =~ s/^\s*//;
        $_ =~ s/\s*$//;
        @g=split(/\s+/,$_);
        $mseek=$g[6];
      }elsif (/WARNING\: Number of orthogonal guesses is/){# considering the H2 case where $mseek is different
        $_ =~ s/^\s*//;
        $_ =~ s/\s*$//;
        @g=split(/\s+/,$_);
        $mseek=$g[7];
      }
    }
  }
  print TSP "\n";
  printf TSF  ("% 18.12f \n"), $Fscurr;
  print NADA "$mseek\n";
  print EN "\n";

  close(EP);
  close(GL);
  close(EN);
  close(NADA);
  close(GD);
  close(FLOG);
  close(GL);

} 

#========================================================================================
#-----------------------------  fromage + Turbomole  ------------------------------------
sub read_turbo_info{
#
#========================================================================================
#
#----------------------------------------------------------------------------------------
#
  my ($methodname,$jobtype);
  %progconf   = prog_config($prog);
  $methodname = $progconf{methodname};

  if ($thres > 0){
    $jobtype="nonadiabatic";
  }

  if (($methodname eq "turbomole-ricc2") or ($methodname eq "turbomole-riadc2")){
    read_cc2_adc2();
  }
  if ($methodname eq "turbomole-tddft"){
    read_dft();
  }
}

#----------------------------------------------------------------------------------------
#
# Read CC2 and ADC2
sub read_cc2_adc2{
  my ($jobtype,$adc,@d1,$is,$ns,$OOS);

  # Read D1 diagnostic
  if ($lvprt >= 3){
    $d1[0]=get_d1("MP2","grad.out");
    if ($d1[0] != -1){print_STDOUT("D1 diagnostic for MP2: $d1[0]\n");}
    $d1[1]=get_d1("CC2","grad.out");
    if ($d1[1] != -1){print_STDOUT("D1 diagnostic for CC2: $d1[1]\n\n");}
    if (($d1[0] > 0.04) or ($d1[1] > 0.05)){
      print_STDOUT("  \n");
      print_STDOUT("   ******************************************\n");
      print_STDOUT("   WARNING: D1 diagnostic indicates that MP2 \n");
      print_STDOUT("   and CC2 may be inadequate to describe the \n");
      print_STDOUT("   ground state for this geometry.           \n");
      print_STDOUT("   ***************************************** \n\n");
    }
  }

  # State information
  print_STDOUT("Information about the excited states:\n",$istep,$kt);
  if ($istep % $kt < $eps){
    my_grep("grad.out","STDOUT"," type: RE0",6,"A","append");
    my_grep("grad.out","STDOUT","contributions of excitation levels to excited states",$nstat+5,"A","append");
  }

  if ($thres > 0){
    $jobtype="nonadiabatic";
  }

  print_STDOUT("Jobtype: $jobtype\n",$istep,$kt);
  # Nonadiabatic coupling
  if ($jobtype eq "nonadiabatic"){
    if ($run_flag ne "second run"){
      if ($lvprt >= 3){printf_STDOUT($istep,$kt,"\nCalling read_write_nonad \n");}
      read_write_nonad();
      if ($lvprt >= 3){printf_STDOUT($istep,$kt,"\n Exiting read_write_nonad \n");}
    }
  }

  # Oscillator strength
  $is = 1;
  for ($ns = 2; $ns <= $nstat; $ns++){
    $OOS = osc_strength($mdle,$is,$ns,$prog);
    if ($lvprt >= 3){printf_STDOUT($istep,$kt,"\nOscillator strength (%d,%d) = %9.6f \n",$is,$ns,$OOS);}
    if (($nxrestart == 0) or ($istep != 0) or ($lvprt >= 2)){
      if ($istep % $kt < $eps){
        printf PP " Oscillator strength (%d,%d) = %9.6f \n",$is,$ns,$OOS;
      }
    }
  }
  print_STDOUT("\n",$istep,$kt);

}
  
#-------------------------------------------------------------------------------------
# Read DFT
sub read_dft{
  my ($typedft,@g,$jobtype,$is,$OOS);

  # State information
  if ($typedft != 3) {
    print_STDOUT("Information about the excited states:\n",$istep,$kt);
    if ($istep % $kt < $eps){
      my_grep("grad.out","STDOUT"," Dominant contributions:",5,"A","append");
    }
  }

  if ($thres > 0){
    $jobtype="nonadiabatic";
  }

  # Nonadiabatic coupling
  if ($jobtype eq "nonadiabatic"){
    if ($run_flag ne "second run"){
      read_write_nonad();
    }
  }

  # Oscillator strength
  $is = 1;
  for ($ns = 2; $ns <= $nstat; $ns++){
    $OOS = osc_strength($mdle,$is,$ns,$prog);
    if ($lvprt >= 3){printf_STDOUT($istep,$kt,"\nOscillator strength (%d,%d) = %9.6f \n",$is,$ns,$OOS);}
    if (($nxrestart == 0) or ($istep != 0) or ($lvprt >= 2)){
      if ($istep % $kt < $eps){
        printf PP " Oscillator strength (%d,%d) = %9.6f \n",$is,$ns,$OOS;
      }
    }
  }
  print_STDOUT("\n",$istep,$kt);

}
#========================================================================================
#========================================================================================
#---------------------------    General Routines ----------------------------------------

sub read_write_nonad{
#
#========================================================================================
#
#----------------------------------------------------------------------------------------
#
  my ($iaux,$jst,$ist,$vh,$kat,$dum);

  # Read and write cioverlap information
  if (( $vdoth == 1 ) or ( $vdoth < 0 )){
    # Read cioverlap information
    if($lvprt >= 2) {print_STDOUT("\nvdoth is $vdoth\n",$istep,$kt);}
    $retval = callprogsp( $mld, "read_cioverlap.pl", $mdle );
    if ($retval != 0){
      die "$mdle is dying now\n";
    }
    # Print cioverlap information
    if ( $lvprt >= 1 ) {print_STDOUT("\nNonadiabatic coupling terms v.h (a.u.):\n",$istep,$kt);}
    open(NV,"nad_vectors") or die "Cannot open nad_vectors.";
    $iaux = 0;
    while(<NV>){
      $iaux++;
      if ( $iaux == 1 ) {
        chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
        ( $vh, $dum, $dum ) = split( /\s+/, $_ );
        printf_STDOUT($istep,$kt,"%14.7F\n", $vh );
        printf_STDOUT($istep,$kt,"%14.7F %14.7F %14.7F\n", $vh, $dum, $dum);
      }elsif( $iaux == $nat ) {
        $iaux = 0;
      }
    }
    close(NV);
  }

} 
