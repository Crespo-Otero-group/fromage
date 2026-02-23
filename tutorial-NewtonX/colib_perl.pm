package colib_perl;
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#====================================================================================
#
#               NEWTON-X PERL LIBRARY COLLECTTION
#               Mario Barbatti 2005 - 2022
#
#====================================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Subroutines and variables exported by this module: 
require Exporter;
@ISA       = qw(Exporter);
@EXPORT    = qw(load_status keyinfile mcscfconv ciconv callprog callprogsp frootci intc changekeyword deletekeyword  
                substitutekeyword  searchkeyword  exit_error  getkeyword  addkeyword  getkeydrt  getkeycolon  number_of_atoms  
                osc_strength  prog_config  rm_files  title_nx  end_nx  max_ns  my_grep  units  print_STDOUT  
                printf_STDOUT  columbus_mode  columbus_mem  test_list  read_control  change_control  delete_control  find_control  get_control  
                print_control  wf_info  make_num_sequence  read_single_value  my_awk save_files get_d1  rewrite_ctd  
                write_internal_control_files  make_gamess_scr read_credits  test_gamess find_scr_gamess load_defaults load_defaults_tpp 
                namelist_leg is_it_numeric addkeyend test_turbomole_version check_thirdparty make_prog_list build_description
		read_cp2kinp change_cp2kinp get_cp2kinp find_cp2kinp find_cp2kinp_full rm_cp2kinp add_cp2kinp add_cp2kkey ciocheck);
@EXPORT_OK = qw($found  $nconv  $niter  $nintcoor  $intcoor  $value  $progname  $method  $description  $dir  $prog  $file  
                $reg_exp  $field  $eps  @vector  $c  $c2  $vec  $more  $nelem  @status  $version);

#
#====================================================================================
#       DATA-INFO ROUTINES
#====================================================================================
#

sub units{
#
#====================================================================================
#
# This subroutine returns the value of a constant, parameter or conversion factor 
# corresponding to a given input variable $key. 
#
#------------------------------------------------------------------------------------
#
    local($key, $value);
    ($key) = @_;
#
    if    ($key eq "BK"      ){$value = 0.3166813639E-5;}   # Boltzmann constant Hartree/kelvin
    elsif ($key eq "bk_ev"   ){$value = 8.617343E-5;}       # Boltzmann constant eV/kelvin
    elsif ($key eq "proton"  ){$value = 1822.888515;}       # (a.m.u.)/(electron mass)
    elsif ($key eq "timeunit"){$value = 24.188843265E-3;}   # femtoseconds
    elsif ($key eq "au2ev"   ){$value = 27.21138386;}       # au to eV
    elsif ($key eq "pi"      ){$value = 3.141592653589793;} # pi
    elsif ($key eq "au2ang"  ){$value = 0.52917720859;}     # au to angstrom
    elsif ($key eq "deg2rad" ){$value = 1.745329252E-2;}    # degree to radian (pi/180)
    elsif ($key eq "au2cm"   ){$value = 219474.625;}        # hartree to cm-1
    elsif ($key eq "cm2au"   ){$value = 4.55633539E-6;}     # cm-1 to hartree
    elsif ($key eq "h_planck"){$value = 1.5198298508E-16;}  # h hartree.s
    elsif ($key eq "h_evs"   ){$value = 4.13566733E-15;}    # h eV.s
    elsif ($key eq "light"   ){$value = 299792458;}         # speed of light m/s
    elsif ($key eq "lightau" ){$value = 137.035999139;}     # speed of light au
    elsif ($key eq "e_charge"){$value = -1.602176487E-19;}  # electron charge C
    elsif ($key eq "e_mass"  ){$value = 9.10938215E-31;}    # electron mass kg
    elsif ($key eq "eps0"    ){$value = 8.854187817e-12;}   # vacuum permitivity C^2*s^2*kg^-1*m^-3
    else{die " Constant $key was not found in colib_perl.pm";}
#
    return $value;
}


sub make_prog_list{
#
#====================================================================================
#
# This subroutine forms program/method description in nxinp.
#
#------------------------------------------------------------------------------------
#
    my ($type,%progconf,$description);
    my ($prog,$ind,$indmax,$key,$label,$method,$menu);

    $progmax = 20;  # If prog goes over 20, increase it here too.

    ($type)=@_; 
    $description = "";

    $indmax = 10*$progmax;
    for ($ind = 0; $ind <= $indmax; $ind++){
         $prog = $ind/10;
         $key    = "";
         %progconf = prog_config($prog);
         $key    = $progconf{key};
         $label  = $progconf{label};
         $method = $progconf{method};
         $menu   = $progconf{$type};
         if ($menu eq "y"){
           if ($key ne ""){
                $description = $description."\t$key - $label";
                if ($method ne "NULL"){
                     $description = $description." $method\n";
                }else{
                     $description = $description."\n";
                }
           }
         }
    } 

    return $description;

}

#
#====================================================================================
#       DYNAMICS AND INITIAL CONDITIONS CONTROLLING
#====================================================================================
#

sub load_status{
#
#====================================================================================
#
# This subroutine returns the values of all control.d variables at a certain time step.
# Usage: 
# ($nat,$istep,$nstat,$nstatdyn,$ndamp,$kt,$dt,$t,$tmax,$nintc, 
#  $mem,$nxrestart,$thres,$killstat,$timekill,$prog,$lvprt,$etot_jump, 
#  $etot_drift,zpecorr)=load_status("control.d",$mdle); 
#            
# Array order:  
#   0:$nat       1:$istep    2:$nstat   3:$nstatdyn   4:$ndamp      5:$kt     6:$dt     
#   7:$t         8:$tmax     9:$nintc  10:$mem       11:$nxrestart  12:$thres 13:$killstat 
#  14:$timekill 15:$prog    16:$lvprt  17:$etot_jump 18:$etot_drift 19:zpecorr
#
# To load particular keyword, for instance $lvprt: 
#  @status = load_status($file,$mdle);
#  $lvprt  = $status[16];
#
#------------------------------------------------------------------------------------
#   
     my ($file,@status,$mdle);
     my ($nat,$istep,$nstat,$nstatdyn,$ndamp,$kt,$dt,$t,$tmax,$nintc,
         $mem,$nxrestart,$thres,$killstat,$timekill,$prog,$lvprt,$etot_jump,
         $etot_drift,$zpecorr);
     $mdle="xx:";
     ($file,$mdle) = @_;
     open(INP,$file) or die "$mdle Cannot open $file!\nsystems answer was: $!\n";
     $_=<INP>;
     chomp;s/\s+//g;
     @status=split(/,/,$_);
     close(INP);
     $istep=$status[1];
     $kt   =$status[5];
     $lvprt=$status[16];
     if ($lvprt>=2) {print_STDOUT("$mdle running STEP=$istep  KT=$kt\n");}
     return @status;
}

sub load_defaults{
#
#====================================================================================
#
# This subroutine contains the defaults for the Newton-X programs. 
#
#------------------------------------------------------------------------------------
#
    my ($param1,$param2);
    my ($kross_d,$cascade_d,$current_d,$never_state_d,$include_pair_d,$e_ci_d,$ci_cons_d,$cio_options_d,$cisc_options_d);
    my ($idalton_d,$cprog_d,$coptda_d,$ncore_d,$ndisc_d,$blasthread_d,$vdoth);
    my ($phase_d,$getphase_d,$ms_d,$integrator_d,$nohop_d,$forcesurf_d,$nrelax_d,$vdoth_d,$tully_d);
    my ($mom_d,$adjmom_d,$probmin_d,$popdev_d,$decohmod_d,$decay_d,$seed_d,$decovlp_d,$iatau_d);
    my ($thrwp_d,$run_complex_d,$gamma_model_d);
    my ($ba_smooth_d,$ba_nsmooth_d,$ba_dh_d,$ba_de_d,$ba_dv_d);
    my ($killstat_d,$timekill_d,$ndamp_d,$lvprt_d,$kt_d,$mem_d,$etot_jump_d,$etot_drift_d);
    my ($thres_d,$nxrestart_d);
    my ($nat_d,$nstat_d,$istep_d,$dt_d,$t_d,$tmax_d,$prog_d);
    my ($numat_d,$nact_d,$iprog_d);
    my ($npoints_d,$fgeom_d,$fout_d,$fvib_d,$chk_e_d,$evert_d,$de_d,$kvert_d,$nis_d,$nfs_d);
    my ($nm_flag_d,$anh_f_d,$rescale_d,$temp_d,$iseed_d,$ics_flg_d,$fnmode_d);
    my ($address_d,$traji_d,$trajf_d,$cmp_e_d,$n_pick_d,$ti_d,$tf_d,$reorder_d,$etot_dev_d,$pop_dev_d);
    my ($ekin_d,$nis_mol_d,$nfs_cat_d,$ip_type_d,$do_type_d,$eq_geom_d,$njobs_dyson_d,$eps_ciscffs_d,$eps_overlap_d);
    my (%progconf,$nad_exec);
    my ($x0_d,$p0_d,$mass_d,$delx_d,$npoints_d);
    my ($subfile_d,$batchdef_d,$title_d,$type_d,$delta_d,$eps_d,$kappa_d,$os_condon_d);
    my ($prob_kind_d,$screen_d,$norm_d,$e_center_d,$e_var_d,$nref_d,$l_shape_d,$dens_d,$absenergy_d,$run_IS_d);
    my ($ics_comp_d,$ezdyson_path_d,$ics_kind_d,$ics_proc_d,$E_points_d,$l_max_d,$Ep_min_d,$Ep_max_d,$Ephoton_d,$Ek_min_d);
    # SM added
    my ($zpecorr_d,$kmodel_d,$kcheck_d,$tavg_d,$tcycle_d,$ethres_d,$biascorrect_d,$kahbond_d,@ind_ah_d,@zperef_d);

    ($file, $param1, $param2) = @_;

    ## Defaults for control.dyn
    #-------------------------------------------------------------------------- 
    if ($file =~ /control/){
       if    ($param1 == 1){
	 $nat_d          = number_of_atoms();
	 $nstat_d        = 2;
	 $istep_d        = 0;
	 $dt_d           = 0.5;
	 $t_d            = 0;
	 $tmax_d         = 100;
	 $prog_d         = sprintf("%.1f",1.0);
	 return $nat_d,$nstat_d,$istep_d,$dt_d,$t_d,$tmax_d,$prog_d;
       }elsif($param1 == 2){
         $killstat_d     = 1;
         $timekill_d     = 0; 
         $ndamp_d        = 0; 
         $kt_d           = 1;      
         $mem_d          = 200;
         $etot_jump_d    = 0.5;   
         $etot_drift_d   = 0.5;
	 $nxrestart_d    = 0;
	 $zpecorr_d      = 0;
	 return $killstat_d,$timekill_d,$ndamp_d,$kt_d,$mem_d,$etot_jump_d,$etot_drift_d,$nxrestart_d,$zpecorr_d;
       }
    }elsif($file =~ /thres/){
       $prog=$param1;
       %progconf   = prog_config($prog);
       $lvprt_d    = $progconf{lvprt_d};
       $nad_exec   = $progconf{nad_exec};
       if ($nad_exec eq "y")  {
	  $thres_d = 100.0;
       }elsif($nad_exec eq "n"){
	  $thres_d = 0.0;
       }
       return $lvprt_d,$thres_d;
    }
    
    ## Defaults for jiri.inp
    #--------------------------------------------------------------------------  
    if ($file =~ /jiri/){
       # Defaults for jiri.inp
       $prog=$param1;
       $vdoth=$param2;
       %progconf   = prog_config($prog);
       $never_state_d  = $progconf{never_state_d};
       $cio_options_d  = $progconf{cio_options_d};
       $cisc_options_d = $progconf{cisc_options_d};
       $cprog_d        = $progconf{cprog_d};
       $kross_d        = 1;
       $cascade_d      = 0;
       $current_d      = 0;
       $include_pair_d = 0;
       $e_ci_d         = 0.2;
       $ci_cons_d      = 1;
       $idalton_d      = 64;
       $coptda_d       = 1;
       $ncore_d        = 0;
       $ndisc_d        = 0;
       $blasthread_d   = 1;
       return $kross_d, $cascade_d, $current_d, $never_state_d, $include_pair_d, $e_ci_d, $ci_cons_d, $cio_options_d,
              $cisc_options_d, $idalton_d , $cprog_d, $coptda_d, $ncore_d, $ndisc_d,$blasthread_d;
    }
    
    ## Defaults for sh.inp
    #--------------------------------------------------------------------------
    if  ($file =~ /sh/){
       $prog     = $param1;
       $vdoth    = $param2;
       %progconf   = prog_config($prog);
       $vdoth_d    = $progconf{vdoth_d};
       $phase_d      = 1;
       if    ($vdoth != 0){
	 if ($progconf{progname} eq "mopac"){      
            $getphase_d   = 1;
         }else{
            $getphase_d   = 0;
	 }
       }elsif($vdoth == 0){
	 if ($progconf{progname} eq "analytical"){
            $getphase_d   = 0;
         }else{
            $getphase_d   = 1;
	 }
       }
       $nohop_d      = 0;
       $forcesurf_d  = 1;
       $nrelax_d     = 0;
       $tully_d      = 1;
       $mom_d        = 1;
       $probmin_d    = 0;
       $decohmod_d   = 1;
       $decay_d      = 0.1;
       $seed_d       = 1;
       $decovlp_d    = 1.0;
       $iatau_d      = 1;
       $thrwp_d      = 0.001;
       if    ($vdoth < 0){
         $ms_d         = 0;
         $integrator_d = 6;
       }elsif($vdoth != -1){
         $ms_d         = 20;
         $integrator_d = 5;
       }
       if    ($vdoth == 0){
         $adjmom_d = 0;
       }elsif($vdoth != 0){
         $adjmom_d = -1;
       }
       $run_complex_d= 0;
       $popdev_d     = 0.05;
       return $ms_d,$integrator_d,$phase_d,$getphase_d,$nohop_d,$forcesurf_d,$nrelax_d,$vdoth_d,$tully_d,$mom_d,$adjmom_d,$probmin_d,$popdev_d,$decohmod_d,$decay_d,$seed_d,$decovlp_d,$iatau_d,$thrwp_d,$run_complex_d;
    }
    
    ## Defaults for auxnac.inp
    #--------------------------------------------------------------------------
    if($file =~ /auxnac/){
       $ba_smooth_d  = 1;
       $ba_nsmooth_d = 4;
       $ba_dh_d      = 0.1;
       $ba_de_d      = 2;
       $ba_dv_d      = 0.1;
       return $ba_smooth_d,$ba_nsmooth_d,$ba_dh_d,$ba_de_d,$ba_dv_d;
    }
    
    ## Defaults for therm.inp
    #--------------------------------------------------------------------------
    if($file =~ /therm/){
       $ktherm_d   = 1; 
       $kts_d      = 1;
       $lts_d      = -1;
       $nstherm_d  = 1;
       $temp_d     = 300;
       $gamma_d    = 0.2;
       $radius_d   = 10;
       $iseed_d    = 1;
       $lvp_d      = 1;
       return $ktherm_d,$kts_d,$lts_d,$nstherm_d,$temp_d,$gamma_d,$radius_d,$iseed_d,$lvp_d;
    }

    ## Defaults for zpe.inp
    #--------------------------------------------------------------------------
    if($file =~ /zpe/){
       $kmodel_d      = 1;
       $kcheck_d      = 0;
       $tavg_d        = 500.0;
       $tcycle_d      = 500.0;
       $ethres_d      = 0.001;
       $biascorrect_d = 0;
       $kahbond_d     = 2;
       @ind_ah_d      = (2,1, 5,4);
       @zperef_d      = (0.0035, 0.0042);
       return $kmodel_d,$kcheck_d,$tavg_d,$tcycle_d,$ethres_d,$biascorrect_d,$kahbond_d,@ind_ah_d,@zperef_d;
    }
    
    ## Defaults for complex.inp (CS-FSSH)
    #--------------------------------------------------------------------------
    if($file =~ /complex/){
       $gamma_model_d     = 1;
       $path_gamma_model_d= "base";
       $prog_gamma_model_d= "my_model";
       $same_mo_d         = 0;
       $ress_shift_d      = 0;
       return $gamma_model_d,$path_gamma_model_d,$prog_gamma_model_d,$same_mo_d,$ress_shift_d;
    }
    
    ## Defaults for initqp_input
    #--------------------------------------------------------------------------
    if($file =~ /initqp_input/){
       if ($param1 == -1){  # get PROG
	    $prog_d    = 6.5;
            return $prog_d;  
       }elsif ($param1 == 0){  # get NACT
            $numat_d   = 3;
            $nact_d    = 2;
            $iprog_d   = 4;
            $ics_flg_d = "n";
            return $numat_d,$nact_d,$iprog_d,$ics_flg_d;  
       }elsif ($param1 == 10){  # get LVPRT
	    $prog      = $param2;
	    if (($prog < 0) or (!defined $prog)){
	      $lvprt_d   = 1;
	    }else{
              %progconf  = prog_config($prog);
              $lvprt_d   = $progconf{lvprt_d};
	    }
            return $lvprt_d;
       }elsif (($param1 == 1) or ($param1 == 2) or ($param1 == 3)){ # NACT=1,2,3
            $npoints_d = 1;
            $fgeom_d   = "geom";
            $fout_d    = "ini_qv";
            $fvib_d    = "qvector";
            $chk_e_d   = 0;
            $evert_d   = 5.0;
            $de_d      = 100.0;
            $kvert_d   = 1;
            $nis_d     = 1;
            $nfs_d     = 2;
            $nm_flag_d = 0;
            $anh_f_d   = 1;
	    if (($param1 == 1) or ($param1 == 2)){
              $rescale_d = "n";
            }elsif($param1 == 3){
              $rescale_d = "y";
	    }
            $temp_d    = 0.0;
            $iseed_d   = -1;
            if ($param2 == 1){$fnmode_d  = "gamess.out";}
            if ($param2 == 2){$fnmode_d  = "force.out";}
            if ($param2 == 3){$fnmode_d  = "suscalls";}
            if ($param2 == 4){$fnmode_d  = "gaussian.log";}
            if ($param2 == 5){$fnmode_d  = "molden.freq";}
            if ($param2 == 6){$fnmode_d  = "freq.out";}
            if ($param2 == 7){$fnmode_d  = "NORMCO";}
            if ($param2 == 8){$fnmode_d  = "";}
            if ($param2 == 9){$fnmode_d  = "modes.xyz";}
            if ($param2 == 11){$fnmode_d  = "orca.hess";}
            return $npoints_d,$fgeom_d,$fout_d,$fvib_d,$chk_e_d,$evert_d,$de_d,$kvert_d,$nis_d,$nfs_d,$nm_flag_d,$anh_f_d,$rescale_d,$temp_d,$iseed_d,$fnmode_d;
       }elsif ($param1 == 4){  # NACT=4
            $address_d = "/home/old_dyn/TRAJECTORIES";
            $traji_d   = 1;
            $trajf_d   = 10;
            $nis_d     = 1;
            $nfs_d     = 2;
            $chk_e_d   = 0;
            $cmp_e_d   = 0;
            $evert_d   = 5.0;
            $de_d      = 0.5;
            $iseed_d   = 0;
            $n_pick_d  = -1;
            $npoints_d  = 1;
            $ti_d       = 0;
            $tf_d       = 0;
            $reorder_d  = 1;
            $etot_dev_d = 0.5;
            $pop_dev_d  = 0.1;
            return $address_d,$traji_d,$trajf_d,$nis_d,$nfs_d,$chk_e_d,$cmp_e_d,$evert_d,$de_d,$iseed_d,$n_pick_d,$npoints_d,$ti_d,$tf_d,$reorder_d,$etot_dev_d,$pop_dev_d;
       }elsif ($param1 == 5){  # NACT=5
            $fgeom_d    = "not found";
            $npoints_d  = 1;
            $ekin_d     = (3*$param2-6)*0.1;
            $temp_d     = 0;
            $iseed_d    = 0;
            return $fgeom_d,$npoints_d,$ekin_d,$temp_d,$iseed_d;
       }elsif ($param1 == 6){  # NACT=6
            $nis_d      = 1;
            $nfs_d      = 2;
            return $nis_d,$nfs_d;
       # SM added    
       }elsif ($param1 == 7){  # NACT=7
            $x0_d      = 0.0;
            $p0_d      = 1.0;
            $mass_d    = 1.0;
            $delx_d    = 0.5;
            $npoints_d = 100;
            $iseed_d   = 0;
            return $x0_d,$p0_d,$mass_d,$delx_d,$npoints_d,$iseed_d;
       }elsif ($param1 == -1){  # ICS_FLG="y"
            $nis_mol_d    = 1;
            $nfs_cat_d    = 1;
            $ip_type_d    = 2;
            $do_type_d    = 1;
            $eq_geom_d    = "n";
            return $nis_mol_d,$nfs_cat_d,$ip_type_d,$do_type_d,$eq_geom_d;
       }
    }
    
    ## Defaults for do.par
    #--------------------------------------------------------------------------
    if($file =~ /do.par/){
       $prog_d = 6.5;
       $njobs_dyson_d = 1;
       $eps_ciscffs_d = 1E-2;
       $eps_overlap_d = 1E-1;
       return $prog_d,$njobs_dyson_d,$eps_ciscffs_d,$eps_overlap_d;
    }

    ## Defaults for mkd.inp
    #--------------------------------------------------------------------------
    if ($file =~ /mkd.inp/){
       if ($param1 eq "Z0"){
         $type_d      = 4;
	 return $type_d;
       }elsif($param1 eq "Z1"){
	 if ($param2 == 2){
           $prob_kind_d = "I";		 
	 }else{	 
           $prob_kind_d = "F";  
         }
	 return $prob_kind_d;
       }elsif($param1 eq "Z2"){	 
         if ($param2 eq "E"){
	   $dens_d      = "Y";
         }else{
	   $dens_d      = "N";
	 }
       }elsif($param1 eq "I"){
         $ics_comp_d    = 1;
         $ezdyson_path  = "/home/bin/ezDyson/exe";
         $ics_kind_d    = 2;
         $ics_proc_d    = 2;
         $E_points_d    = 2;
         $l_max_d       = 6;
         $Ep_min_d      = 5.0;
         $Ep_max_d      = 20.0;
         $Ephoton_d     = 10.0;
         $Ek_min_d      = 0.0;	      
	 return $ics_comp_d,$ezdyson_path_d,$ics_kind_d,$ics_proc_d,$E_points_d,$l_max_d,$Ep_min_d,$Ep_max_d,$Ephoton_d,$Ek_min_d; 
       }else{
         $subfile_d   = "pmold";
         $batchdef_d  = "#PBS -N";
         $title_d     = "title";
         $seed_d      = 0;
         $nis_d       = 1;
         $nfs_d       = 2;
         $delta_d     = 0.01;
         $eps_d       = 0.005;
         $kappa_d     = 0;
         $os_condon_d = -1;
         $screen_d    = 0;
         $norm_d      = "local";
         $e_center_d  = 0;
         $e_var_d     = 0.5;
         $nref_d      = 1;
         $temp_d      = 0;
         $l_shape_d   = "gauss";
         $absenergy_d = 0;
         $run_IS_d    = 0;
         return $subfile_d,$batchdef_d,$title_d,$seed_d,$nis_d,$nfs_d,$delta_d,$eps_d,$kappa_d,$os_condon_d,$screen_d,$norm_d,$e_center_d,$e_var_d,$nref_d,$temp_d,$l_shape_d,$absenergy_d,$run_IS_d;
       }
    }

}

#
#====================================================================================
#       PROGRAM-EXECUTION ROUTINES
#====================================================================================
#

sub callprog{
#
#====================================================================================
#
# This subroutine runs executable $prog (with path to it $mld) and checks for 
# successful execution. It does basically the same as 'callprogsp'.      
#    
# Usage:   callprog($mld, $prog, $mdle, $BASEDIR)         
# - If $mld is included in the environment variable $PATH, enter $mld = "".     
# - $mdle is a flag, usually with the program name from which callprogsp is called.   
# - $BASEDIR: runtime error will be written into $BASEDIR/../DEBUG/runnx.error file.
#
#------------------------------------------------------------------------------------
#
    local ($mld, $mdle, $x, $prog, $scrv, $status, $t, $istep);
    $mld     = $_[0];
    $prog    = $_[1];
    $mdle    = $_[2];
    $BASEDIR = $_[3];
    if (! $BASEDIR){
       $BASEDIR = `pwd`;
       chomp ($BASEDIR);
    }
    $x = $prog;
    $x =~ s/[\<\>].*$//;
#
    print_STDOUT("---------------------------------------------------------------------------\n");
    $started = "Starting $x at " . qx/date/;
    print_STDOUT("$started");
    if ($mld eq ""){
       $scrv = system("$prog 2>> $BASEDIR/../DEBUG/runnx.error");
    }else{
       $scrv = system("$mld/$prog 2>> $BASEDIR/../DEBUG/runnx.error");
    }
    $status = "";
    if ($scrv == 0){
       if ($x =~ m/runc/){
          lookatcolumbus("runc.error");
          lookatcolumbus("runls");
       }elsif ($x =~ m/gau/){
          lookatgau("gaussian.log");
       }
    }else{
       if($x =~ m/gau/){
          lookatgau("gaussian.log");
       }else{
          $status = "with ERROR";
          if ($x =~ m/runc/){
             test_columbus_version();
          }
          err_message();
       }
    }
    if ($status ne "with ERROR"){
       $status = "successfuly";
    }
    $finished = "Finished $x $status at " . qx/date/;
    print_STDOUT("$finished");
    print_STDOUT("---------------------------------------------------------------------------\n");
    $! = 0;
    if ($status eq "with ERROR"){
       $! = 256;                    # for the moment error number 256 is used for anything
       if ($mdle eq "moldyn.pl:") {
          cptemp();
          print_STDOUT("See also error messages in DEBUG/runnx.error \n");
       }
    }
    return $!;
}

sub callprogsp{
#
#====================================================================================
#
# This subroutine runs executable $prog (with path to it $mld) and checks for 
# successful execution. Call it as: 
#                                                                                                    
# Usage:   callprog($mld, $prog, $mdle, $BASESIR)                          
# - If $mld is included in the environment variable $PATH, enter $mld = "". 
# - $mdle is a flag, usually with the program name from which callprogsp is called. 
# - $BASEDIR: runtime error will be written into $BASEDIR/../DEBUG/runnx.error file.
#
#------------------------------------------------------------------------------------
#
    local ($mld, $mdle, $x, $prog, $scrv, $status, $t, $istep);
    $mld     = $_[0];
    $prog    = $_[1];
    $mdle    = $_[2];
    $BASEDIR = $_[3];
    if (!$BASEDIR){
       $BASEDIR = `pwd`;
       chomp($BASEDIR);
    }
    $x = $prog;
    $x =~ s/[\<\>].*$//;
    $started = "Starting $x at " . qx/date/;
# 
    if ($mld eq ""){
       $scrv = system("$prog 2>> $BASEDIR/../DEBUG/runnx.error");
    }else{
       $scrv = system("$mld/$prog 2>> $BASEDIR/../DEBUG/runnx.error");
    }
    $status = "";
    if ($scrv == 0){
       if ($x eq "runc"){
          lookatcolumbus("runc.error");
          lookatcolumbus("runls");
       }
       if ($status ne "with ERROR"){
          $status = "successfuly";
       }
    }else{
       $status = "with ERROR";
       err_message();
    }
#
    $finished = "Finished $x $status at " . qx/date/;
    $! = 0;
    if ($status eq "with ERROR"){
       $! = 256;                    # for the moment error number 256 is used for anything
       if ($mdle eq "moldyn.pl:"){
          cptemp();
          print_STDOUT("---------------------------------------------------------------------------\n");
          print "$started";
          print_STDOUT("$finished");
          print_STDOUT("See also error messages in DEBUG/runnx.error \n");
          print_STDOUT("---------------------------------------------------------------------------\n");
       }
    }
    return $!;
}

#
#====================================================================================
#       REGEX AND TEST ROUTINES
#====================================================================================
#

sub keyinfile{
#
#====================================================================================
#
# Given an namelist-type input file ($filename), i.e., each line containing 
# $key1 = valor1[,], $key2 = valor2[,], ... and a variable ($key), this subroutine 
# searches $key in $filename and returns its value ($found) if found. If not, returns 
# the empty value $found = "".
#
#------------------------------------------------------------------------------------
#
    local ($i, $j, $k, @keysperline, $found, @keyval, $filename, $key);
    ($filename, $key) = @_;
    local ($finished);
#
    open (ANYFILE, "$filename");
    $/ = "\n";
    @keyval = <ANYFILE>;
    close ANYFILE;
    $found = "";
    foreach $i (@keyval){
         chop $i;
         $i =~ s/ *$//g;
         @keysperline = split ',',$i;
         $k = 0;
         $finished = 0;
         for ($j = 0; $j <= $#keysperline; $j++){  
             $_ = $keysperline[$j];
             if (/\b$key\b/i) {$found = $keysperline[$j]; $k= 1;next;}
             if ( ( ! /[\$&A-Z\/]/i ) && $k ) { $found = join ':', $found, $keysperline[$j];}
             if (( /[\$&A-Z\/]/i ) && $k ) {$finished = 1;}
         }
         if ($finished) {last;}
    }
    $found =~ s/^.*= *//g;
    chomp($found);
    $found =~ s/\'//g;
    $found =~ s/\"//g; # clean quotation marks
    return $found;
}

sub changekeyword{
#
#====================================================================================
#
# Given two namelist-type input files ($file_old and $file_new), a keyword ($keyword) 
# and an input value ($value) for $keyword, this subroutine copies the content of 
# $file_old into $file_new, with the corresponding value of $keyword replaced by $value. 
# If $file_old and $file_new are the same file, simply replaces the value of $keyword 
# by $value.
#
# Usage:                                           
# changekeyword($file_old, $file_new, $keyword, $value)    
# 
#------------------------------------------------------------------------------------
#
    local ($file_old, $file_new, $keyword, $value, $found);
    ($file_old, $file_new, $keyword, $value) = @_;
#
    $/ = "\n";
    open OLD, "<$file_old";
    open TMP, ">tmp-clp";
    while (<OLD>){print TMP $_;}
    close OLD; 
    close TMP;
#
    open OLD, "<tmp-clp";
    open NEW, ">$file_new";
    $found = 0;
    while (<OLD>){
	  if (/\b$keyword\b/i){ 
             print NEW "  $keyword = $value\n";
             $found = 1; 
          }else{ 
             print NEW $_;
          }
    }
    close NEW;
    close OLD;
    if (! $found == 1){die "keyword $keyword not found in file: $file_old\n";}
    unlink("tmp-clp");
    return 0;
}

sub deletekeyword{
#
#====================================================================================
#
# Delete keyword
#
#------------------------------------------------------------------------------------
#
       local ($file_old,$file_new,$keyword,$found);
       ($file_old,$file_new,$keyword)=@_;
#
       $/="\n";
       open OLD, "<$file_old";
       open TMP, ">tmp-clp";
       while(<OLD>) {print TMP $_;}
       close OLD; close TMP;
#
       open OLD, "<tmp-clp";
       open NEW, ">$file_new";
       $found=0;
       while ( <OLD> ){
         if (/\b$keyword\b/i) { 
            $found=1; 
         }else{ 
            print NEW $_;
         }
       }
       close NEW;
       close OLD;
       if (! $found == 1){ die "keyword $keyword not found in file: $file_old\n";}
       return 0;
      }

sub substitutekeyword{
#
#====================================================================================
#
# Substiture keyword.
#
#------------------------------------------------------------------------------------
#
       local ($file_old,$file_new,$keyword,$keyword_new,$value,$value_new,$found);
       ($file_old,$file_new,$keyword,$keyword_new,$value_new)=@_;
#
       $/="\n";
       open OLD, "<$file_old";
       open TMP, ">tmp-clp";
       while(<OLD>) {print TMP $_;}
       close OLD; close TMP;
#
       open OLD, "<tmp-clp";
       open NEW, ">$file_new";
       $found=0;
       while ( <OLD> )
        {
         if (/\b$keyword\b/i) { print NEW "  $keyword_new=$value_new\n";$found=1; }
         else         { print NEW $_;}
        }
       close NEW;
       close OLD;
       if (! $found == 1){ die "keyword $keyword in not found in file: $file_old\n";}
       return 0;
      }

sub searchkeyword{
#
#====================================================================================
#
# Given an input file ($filename) and a keyword ($keyword), this subroutine searches 
# $keyword in $filename and returns the number of lines $keyword is present in the 
# file ($found), whether with upper or lower case.    
#  
# Usage:        
# $result = searchkeyword($filename, $keyword) 
#
#------------------------------------------------------------------------------------
#
    local ($file, $keyword, $found);
    ($file, $keyword) = @_;
    $found = 0;
#
    if (-s $file){
       open(FILE, "<$file") or warn "Cannot read $file !";
       while(<FILE>){
          chomp $_;
          if (/\b$keyword\b/i){$found = $found + 1;}
       }
       close FILE;
    }
    $keyword = $found;
    return $found;
}

sub getkeyword{
#
#====================================================================================
#
# Given a namelist-type input file ($file) and a keyword ($keyword), this subroutine 
# searches $keyword in $filename and returns its value ($value) if found, otherwise 
#returns an input value stored in $default ($value = $default). 
#
# Usage: 
# $result = getkeyword($file, $keyword, $default)   
#
#------------------------------------------------------------------------------------
#
    local ($file, $keyword, $default, $found, $value);
    ($file, $keyword, $default) = @_;
    $found = -1;
    if (-s $file){
       $found = 0;
       $found = searchkeyword("$file", "$keyword");
       if ($found == 0){
          $value = $default;
       }else{
          $value = keyinfile("$file", "$keyword");
       }
    }else{
       $value = $default;
    }
    return $value;
}

sub getkeycolon{
#
#====================================================================================
#   
# Search for the first occurence of $key in $file.
# Return $value as the remaing after colon.
# If not found, return $default
# Example:
# If LISTINGS/cidrtls.all contains the line
# total molecular orbitals            :    36
# then nmol = 36 when the subroutine is invoked in this way
# $nmol= getkeycolon("LISTINGS/cidrtls.all","total molecular orbital","0");
#
#------------------------------------------------------------------------------------
#
 my ($file,$key,$default,$value,$dum);
 ($file,$key,$default)=@_;
 open(FL,$file) or die ":( $file";
 while(<FL>){
   $value=$default;
   if (/\b$key\b/i){
      chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
      ($dum,$value)=split(/:/,$_);
       $value =~ s/^\s*//;$value =~ s/\s*$//;
       last;
   }
 }
 close(FL);
 return $value;
}

sub addkeyword{
#
#====================================================================================
#
# This subroutine adds a keyline (key = value) to a file on the first possible 
# position. Use it as:
#     
# addkeyword(old_filename, new_filename, key, value)                                                           
#                                                                                                               
# - if old_filename = new_filename  ->  the file is overwritten    
# - if value is omitted             ->  only the keyword is written     
# - if the key is already present   
#   ->  either nothing is done (if it is a single key)
#   ->  changekeyword($file_old, $file_new, $keyword, $value) is used to change value
#   ->  a single key will replace a key = value pair     
#   ->  a key = value pair will replace a single key  
#
#------------------------------------------------------------------------------------
#
    local($file_old, $file_new, $keyword, $value, $found);
    my $mdle = "addkeyword" ;
# Got three arguments indicating that a line with $keyword has to be inserted (this cannot be a namelist file). 
# - $keyword is added in first line for empty files and in second line of files with contents
    if ($#_ == 2){
       ($file_old, $file_new, $keyword) = @_;
       saveoldfile($file_old);
       $found = searchkeyword("tmp-clp", $keyword);
       if ($found == 0){
# ..Key not existent in old file - adding it in new file
          open TMP, "<tmp-clp" or open_err("tmp-clp", "read", "die", $mdle);
          open NEW, ">$file_new" or open_err($file_new);
          if (-z "tmp-clp"){
             print NEW "$keyword";
          }else{
             $_ = <TMP>;
             print NEW $_;
             print NEW "$keyword\n";
             while(<TMP>){
                  print NEW $_;
             }
          }
          close TMP or close_err("tmp-clp");
          close NEW or close_err($file_new);
       }else{
# ..Key exists in old file - copy new file
          open TMP, "<tmp-clp" or open_err("tmp-clp");
          open NEW, ">$file_new" or open_err($file_new);
          while (<TMP>){
# ...if-else neccesary to substitute key = value by key only
                if (/\b$keyword\b/){
                   print NEW "$keyword\n";
                }else{
                   print NEW $_;
                }
          }
          close TMP or close_err("tmp-clp");
          close NEW or close_err($file_new);
       }
    }
# Got four arguments indicating thet a line with "key = value" has to be inserted or the value changed.
# (this can be a namelist file or not)
# - $keyword = value is added in first line for empty files and in second line to files with contents
    if ($#_ == 3){
       ($file_old, $file_new, $keyword, $value) = @_;
       saveoldfile($file_old);
       $found = searchkeyword("tmp-clp", $keyword);
       if ($found == 0){
# ..key not existent in old file - adding it in new file
          open TMP, "<tmp-clp" or open_err("tmp-clp");
          open NEW, ">$file_new" or open_err($file_new);
          if (-z "tmp-clp"){
             print NEW " $keyword = $value\n";
          }else{
             $_ = <TMP>;
             print NEW $_;
             print NEW " $keyword = $value\n";
             while (<TMP>){
                   print NEW $_;
             }
          }
          close TMP or close_err("tmp-clp");
          close NEW or close_err($file_new);
       }else{
# ..key exists in old file - use changekeyword()
          changekeyword($file_old,$file_new,$keyword,$value);
       }
    }
    system("rm -f tmp-clp");
}

sub addkeyend{
#
#====================================================================================
#
# This subroutine adds a keyline (key = value) to a file at the last possible 
# position. Use it as:   
# addkeyend(old_filename, new_filename, key, value)
# 
# - if old_filename = new_filename  ->  the file is overwritten 
# - if the key is already present   
#   ->  changekeyword($file_old, $file_new, $keyword, $value) is used to change value
#   ->  a single key will replace a key = value pair
#   ->  a key = value pair will replace a single key
#
#------------------------------------------------------------------------------------
#
  local($file_old, $file_new, $keyword, $value, $found, $end_mark, $line);
  my $mdle = "addkeyend";
  ($file_old, $file_new, $keyword, $value) = @_;
  saveoldfile($file_old);
  $found = searchkeyword("tmp-clp", $keyword);
  $end_mark="n";
  if ($found == 0){
# ..key not existent in old file - adding it in new file
     open TMP, "<tmp-clp" or open_err("tmp-clp");
     open NEW, ">$file_new" or open_err($file_new);
     if (-z "tmp-clp"){
        print NEW " $keyword = $value\n";
     }else{
        while (<TMP>){
           chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
           $line=$_;
           if ($line=~/^\//){
             $end_mark = "y";
           }elsif($line=~/&end/i){
             $end_mark = "y";
           }elsif($line eq ""){
             $end_mark = "y";
           }elsif(!defined $line){
             $end_mark = "y";
           }elsif(eof(TMP)){
             $end_mark = "eof";
           }
           if ($end_mark eq "n"){
              if (($line=~/^&/) or ($line=~/^\//)){
                print NEW "$line\n";
              }else{
                print NEW " $line\n";
              }
           }elsif($end_mark eq "y"){
              print NEW " $keyword = $value\n";
              $end_mark = "n";
              if (($line=~/^&/) or ($line=~/^\//)){
                print NEW "$line\n";
              }else{
                print NEW " $line\n";
              }
           }elsif($end_mark eq "eof"){
              if (($line=~/^&/) or ($line=~/^\//)){
                print NEW "$line\n";
              }else{
                print NEW " $line\n";
              }
              print NEW " $keyword = $value\n";
              $end_mark = "n";
           }
        }
     }
     close TMP or close_err("tmp-clp");
     close NEW or close_err($file_new);
  }else{
# ..key exists in old file - use changekeyword()
     changekeyword($file_old,$file_new,$keyword,$value);
  }
  system("rm -f tmp-clp");
}

sub read_single_value{
#
#====================================================================================
#
# Reads the first space separated value in a file and return it.
# Usage:
# $value=read_single_value(file,default);
#
#------------------------------------------------------------------------------------
#
  local ($value,$file,$default);
  ($file,$default)=@_;
  if (-s $file){
    open(FL,$file) or die ":( $file\n";
    $_=<FL>;
    chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
    $value=$_;
    close(FL);
  }else{
    $value = $default;
  }
  return $value;
}

sub my_grep{
#
#====================================================================================
#
# This routine works like the system grep with options -A and -B.
# Usage:
# my_grep($file,$file_dest,$pattern,$lines,$direction,$position)
# where $file - file to be searched
#       $file_dest - destination file
#                    "STDOUT" = print to standard output
#       $pattern - the pattern to be searched
#       $lines - number of lines to be printed before or after the matching (0 = default)
#       $direction A - print $lines after matching (default)
#                  B - print $lines before matching
#       $position new - overwrite previous content (default)
#                 append - append to previous content
#
#------------------------------------------------------------------------------------
# 
 local($file,$file_dest,$pattern,$lines,$direction,$position,$found,$i,$j,@history,$mdle);
 $mdle  = "my_grep subroutine: ";
 $lines = 0;
 $direction = "A";
 $position = "new";
 $found = "";
 $file = "none_in_mygrep";
 $file_dest = "none_in_mygrep";
 ($file,$file_dest,$pattern,$lines,$direction,$position)=@_;
  $direction = lc $direction;
  $position  = lc $position;
  open(FL,$file) or die "Cannot open $file!";
  while(<FL>){
    push(@history,$_);
    $j++;
    if (/$pattern/){
      $i=0;
      $found = $found.$_;
      if ($lines != 0){
        if ($direction eq "a"){
          while(<FL>){
            $i++;
            if ($i <= $lines){
              $found = $found.$_;
            }elsif($i == $lines+1){
              $found = $found."\n";
              last;
            }
          }
        }elsif($direction eq "b"){ # Not tested!
          for ($i = 1; $i<=$lines; $i++){
            $found = $history[$j-$i-1].$found;
            if ($j-$i-1 <= 0){
              last;
            }
          }
        }
      }
    }
  }
  close(FL);
  if ($file_dest eq "STDOUT"){
    print_STDOUT($found);
  }else{
    if ($position eq "new"){
      open(FL,">$file_dest") or die "$mdle Cannot open $file_dest to write!";
    }elsif($position eq "append"){
      open(FL,">>$file_dest") or die "$mdle Cannot open $file_dest to append!";
    }else{
      die "$mdle exiting without writing \n";
    }
    print FL $found;
    close(FL);
  }
}

sub my_awk{
#
#====================================================================================
#
# This routine works like awk '/reg_exp/ {print $(field+1)}' file
# This means, it looks for a regular expression and returns the argument numebr "field".
# (First argument is 0.)
# If position = "first", the first occurence is returned.
# if position = "last", last occurence is returned.
# If reg_exp is not found, "NOT FOUND" is returned.
#
# Use:
# $value=my_awk($file,$reg_exp,$field,$position);
#
#------------------------------------------------------------------------------------
#
  my ($file,$reg_exp,$field,$position,$value,@g);
  ($file,$reg_exp,$field,$position)=@_;
  $value = "NOT FOUND";
  if (-s $file){
    open(FL,$file);
    while(<FL>){
      if (/$reg_exp/){
         chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
         @g=split(/\s+/,$_);
         if ($position eq "first"){
           last;
         }
      }
    }
    close(FL);
    $value=$g[$field];
  }
  return $value;
}

sub namelist_leg{
#
#====================================================================================
#
# This subroutine rewrites namelists in standard formats.  
# 
# Usage:   
# namelist_leg(<file name>,<print level>);
#
#------------------------------------------------------------------------------------
#
  my ($inpf,$lvprt,$typ,@g);
  ($inpf,$lvprt)=@_;
  if ($lvprt >= 2){
    print "namelist_leg: $inpf will be formated to comply to fortran standards. \n";
  }
  open(JIN,"$inpf");
  open(JINA,">nmlleg-aux");
  while(<JIN>){
    chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
    if ($_ eq ""){                                  # if empty, do nothing    
    }elsif(/&/ and !/end/i){                        # print head line
      print JINA "$_\n";
    }elsif(!/&/ and !/^\//){                         # print variable lines
      @g=split(/=/,$_);
      $g[0] =~ s/^\s*//;$g[0] =~ s/\s*$//;
      $g[1] =~ s/^\s*//;$g[1] =~ s/\s*$//;
      $typ=is_it_numeric($g[1]);
      if ($typ == 0){                               # numerical variables
        print JINA " $g[0] = $g[1]\n";
      }elsif($typ == 1){                            # non-numerical variables
        if (/"/){                                   # with quotation
           print JINA " $g[0] = $g[1]\n"; 
        }else{                                      # add quotation
           print JINA " $g[0] = \"$g[1]\"\n"; 
        }
      }
    }elsif(/&end/ or /&END/){                       # print closing: if &end, replace by /
      print JINA "/ \n";
    }elsif(/\//){                                   # print default closing /
      print JINA "$_ \n";                          
    }
  }
  close(JINA);
  close(JIN);
  system("cp -f nmlleg-aux $inpf; rm -f nmlleg-aux");
}

sub intc{
#
#====================================================================================
#
# Read intcfl
#
#------------------------------------------------------------------------------------
#
      my (@coortype,$i,@field,$string);
      @coortype=qw(stre bend out tors lin1 lin2 tor2);
      @hessdiag=();
      $/="\n";
      $intcoor="";
      open(INTC,"intcfl") or die(" File: intcfl missing; generate internal coordinates first!");
      $title=<INTC>;
      $i=1;
      $nintcoor=0;
      $icount=0;
      while (<INTC>) {
      $line=$_;
      $line=~tr/A-Z/a-z/;
      @field=split('');
      $string=join'',@field[20..23];
      $kstring=$field[0];
      $string=~tr/A-Z/a-z/;
      $string=~s/\W.*//;
      $kstring=~s/\W.*//;
      $string=~tr/A-Z/a-z/;
      $kstring=~tr/A-Z/a-z/;
        for ($i=0; $i<=7;$i++)
         {
          if ($i==7)
           {
            # read the force constat diagonals
            $line=~s/^ *//;
            (@tmp)=split(/\s+/,$line,8);
            @hessdiag= (@hessdiag,@tmp);
           }
          if ($string eq $coortype[$i] && $kstring eq "k"){$intcoor=$intcoor.$line;$nintcoor++;last}
          if ($string eq $coortype[$i]){$intcoor=$intcoor.$line;last}
         } # end of for ($i=0; $i<=6;$i++)
      } # end of: while (<INTC>)
      close(INTC);
     return $nintcoor,$intcoor,\@hessdiag ;
  } 

sub number_of_atoms{
#
#====================================================================================
#
# This subroutine counts the number of atoms according to the information given in 
# at least one of the files: geom, initial_condition or ../geom and returns this value 
# as $value. If non of the files exist, they are empty or cannot be opened, returns $value = 0.
#
#------------------------------------------------------------------------------------
#
    local ($file, $file2, $value);
    $file = "geom";
    if (-s "initial_condition"){
       $file2 = "initial_condition";
    }else{
       $file2 = "initial_condition.old";
    }
    $file3 = "../geom";
    $value = 0;
    if (!-s $file){
       if (-s $file2){
          open(FL, $file2) or return $value;
          FL_LINES: while (<FL>){
                if (/geometry in/i){
                   $value = 0;
                   while (<FL>){
                         $value++;
                         if (/velocity/i){
                            $value = $value - 1;
                            last FL_LINES;
                         }
                   }
                }
          }
          close(FL);
       }else{
          if (-s $file3){
             open(FL, $file3) or return $value;
             while (<FL>){
                   $value++;
             }
             close(FL);
             return $value;
          }
       }
       return $value;
    }else{
       open(FL,$file) or return $value;
       while (<FL>){
             $value++;
       }
       close(FL);
       return $value;
    }
}

sub max_ns{
#
#====================================================================================
#
# This subroutine determines the maximum value between 'nfs' (final electronic state) 
# and 'nis' (initial electronic state).
#
#------------------------------------------------------------------------------------
#
    local ($nis, $nfs, $value);
    ($nis, $nfs) = @_;
    if ($nis > $nfs){
       $value = $nis;
    }elsif($nis < $nfs){
       $value = $nfs;
    }elsif($nis == $nfs){
       $value = 0;
       exit_error("NIS cannot be equal NFS. Check input files.".$mdle);
    }
    return $value;
}

sub make_num_sequence{
#
#====================================================================================
#
# This subroutine takes a positive integer sentence like "1,5-7,10-13,9,15" and 
# returns: (1,5,6,7,9,10,11,12,13,15). Returns error |
# if sentence contains A-Z characters. Redundancies are eliminated. Initial sentence 
# do not need to be sorted.                     |
#
# Usage:  $line="1,5-7,10-13,9,15";    
#         @my_sequence = make_num_sequence($line);  
#
#------------------------------------------------------------------------------------
#
    my ($line, @vector, $i, @a, $k, @s, $last_seen);
    ($line) = @_;
    chomp($line);
    $line =~ s/\s+//g;                          # elimilate spaces
    if ($line =~ /[A-Z]/){
       die "ND\n";                              # check whether it contain A-Z charac.
    }
    @vector = split(/,/, $line);                # split at comma
    $i = 0;
    foreach(@vector){
       if (/-/){
          @s = split(/-/, $_);                  # split at dash
          for ($k = $s[0]; $k <= $s[1]; $k++){
              $a[$i] = $k;                      # accumulates sequence
              $i++;
          }
       }else{
          $a[$i] = $_;                          # accumulates sequence
          $i++;
       }
    }
    @a = sort{$a<=>$b} @a;                      # sort numerical array
    $last_seen = -1.1;
    $i = 0;
    foreach(@a){                                # eliminate redundancies
       if ($_ != $last_seen){
          $vector[$i] = $_;
          $i++;
       }
       $last_seen = $_;
    }
    return @vector;
}

sub rewrite_ctd{
#
#====================================================================================
#
# Do the task of "inp". Use getkeyword to get the values from control.dyn
# and set the defaults at the same time.
# New keywords do not break everything any more.
# Writing control.d and control.d1 is in a separate routine which is
# also called from the hybrid_prep_*.pl.
#
#------------------------------------------------------------------------------------
#
  my ($nat, $istep, $nstat, $nstatdyn, $ndamp, $kt, $nxrestart, $nintc, $mem);
  my ($killstat, $lvprt, $dt, $t, $tmax, $thres, $timekill, $prog, $etot_jump, $etot_drift, $zpecorr);
  my ($line, @write_pos, @write_vel);
  my ($nat_d,$nstat_d,$istep_d,$dt_d,$t_d,$tmax_d,$prog_d);
  my ($killstat_d,$timekill_d,$ndamp_d,$lvprt_d,$kt_d,$mem_d,$etot_jump_d,$etot_drift_d,$nxrestart_d,$zpecorr_d);
  my ($thres_d);
  my $ctdyn = "control.dyn";

   # Defaults for dynamics control (1)
  ($nat_d,$nstat_d,$istep_d,$dt_d,$t_d,$tmax_d,$prog_d)=load_defaults($ctdyn,1,"");

  $nat        = getkeyword($ctdyn,"nat",$nat_d);
  $nstat      = getkeyword($ctdyn,"nstat",$nstat_d);
  $nstatdyn   = getkeyword($ctdyn,"nstatdyn",$nstat);
  $istep      = getkeyword($ctdyn,"istep",$istep_d);
  $dt         = getkeyword($ctdyn,"dt",$dt_d);
  $t          = getkeyword($ctdyn,"t",$t_d);
  $tmax       = getkeyword($ctdyn,"tmax",$tmax_d);
  $prog       = getkeyword($ctdyn,"prog",$prog_d);

  # Defaults for dynamics control (2)
  ($killstat_d,$timekill_d,$ndamp_d,$kt_d,$mem_d,$etot_jump_d,$etot_drift_d,$nxrestart_d,$zpecorr_d)=load_defaults($ctdyn,2,"");
  ($lvprt_d,$thres_d)=load_defaults("thres",$prog,"");

  $thres      = getkeyword($ctdyn,"thres",$thres_d);
  $killstat   = getkeyword($ctdyn,"killstat",$killstat_d);
  $timekill   = getkeyword($ctdyn,"timekill",$timekill_d);
  $ndamp      = getkeyword($ctdyn,"ndamp",$ndamp_d);
  $lvprt      = getkeyword($ctdyn,"lvprt",$lvprt_d);
  $kt         = getkeyword($ctdyn,"kt",$kt_d);
  $mem        = getkeyword($ctdyn,"mem",$mem_d);
  $etot_jump  = getkeyword($ctdyn,"etot_jump",$etot_jump_d);
  $etot_drift = getkeyword($ctdyn,"etot_drift",$etot_drift_d);
  $nxrestart  = getkeyword($ctdyn,"nxrestart",$nxrestart_d);
  $zpecorr    = getkeyword($ctdyn,"zpecorr",$zpecorr_d);

  $nintc     = 3*$nat-6;
  &write_internal_control_files($nat, $istep, $nstat, $nstatdyn, $ndamp, $kt, $dt, $t, $tmax, $nintc, $mem, $nxrestart, $thres, $killstat, $timekill, $prog, $lvprt, $etot_jump, $etot_drift, $zpecorr);

  # Write new control.dyn
  open(IN,">$ctdyn") or die "Cannot write to $ctdyn";
  print IN " &input\n";
  print IN "\tnat        = $nat\n";
  print IN "\tnstat      = $nstat\n";
  print IN "\tnstatdyn   = $nstatdyn\n";
  print IN "\tistep      = $istep\n";
  print IN "\tt          = $t\n";
  print IN "\tdt         = $dt\n";
  print IN "\ttmax       = $tmax\n";
  print IN "\tprog       = $prog\n";
  print IN "\tthres      = $thres\n";
  print IN "\tkillstat   = $killstat\n";
  print IN "\ttimekill   = $timekill\n";
  print IN "\tndamp      = $ndamp\n";
  print IN "\tlvprt      = $lvprt\n";
  print IN "\tkt         = $kt\n";
  print IN "\tmem        = $mem\n";
  print IN "\tetot_jump  = $etot_jump\n";
  print IN "\tetot_drift = $etot_drift\n";
  print IN "\tnxrestart  = $nxrestart\n";  
  print IN "\tzpecorr    = $zpecorr\n";
  print IN "/\n";  
  close(IN);
}

sub write_internal_control_files{
#
#====================================================================================
#
# write the formatted control.d and control.d1
#
#------------------------------------------------------------------------------------
#
  my ($nat, $istep, $nstat, $nstatdyn, $ndamp, $kt, $dt, $t, $tmax, $nintc, $mem, $nxrestart, $thres, $killstat, $timekill, $prog, $lvprt, $etot_jump, $etot_drift, $zpecorr)=@_;

  open(CTD,">control.d") or die "cannot write control.d\n$!\n";
  open(CTD1,">control.d1") or die "cannot write control.d1\n$!\n";

  printf CTD "%8d,%8d,%8d,%8d,%8d,%8d,%12.4f,%12.4f,%12.4f,%8d,%8d,%8d,%8.2f,%8d,%12.4f,%9.1f,%8d,%8.3f,%8.3f,%8d\n", 
  $nat, $istep, $nstat, $nstatdyn, $ndamp, $kt, $dt, $t, $tmax, $nintc, $mem, $nxrestart, 
  $thres, $killstat, $timekill, $prog, $lvprt, $etot_jump, $etot_drift, $zpecorr;

  printf CTD1 "%8d  %8d  %8d  %8d  %8d  %8d  %12.4f  %12.4f  %12.4f  %8d  %8d  %8d  %8.2f  %8d  %12.4f  %9.1f  %8d  %8.3f  %8.3f  %8d\n", 
  $nat, $istep, $nstat, $nstatdyn, $ndamp, $kt, $dt, $t, $tmax, $nintc, $mem, $nxrestart, 
  $thres, $killstat, $timekill, $prog, $lvprt, $etot_jump, $etot_drift, $zpecorr;

  close(CTD) or warn "could not close control.d properly\n$!\n";
  close(CTD1) or warn "could not close control.d1 properly\n$!\n";
}
 
sub is_it_numeric{
#
#====================================================================================
#
# This subroutine returns 0 if input is numeric and 1 if it is not. 
#
# Usage:
# $label=is_it_numeric(<input>);
#
#------------------------------------------------------------------------------------
#
  my ($input,$value);
  # Is input numeric?
  ($input)=@_;
  if ( $input =~ /^[\+-]*[0-9]*\.*[0-9]*$/ && $input !~ /^[\. ]*$/  ) {
    $value=0;  #numeric
  }else{
    $value=1;  #non-numeric
  };
  return $value
}

#
#====================================================================================
#       OUTPUT-CONTROL ROUTINES
#====================================================================================
#

sub title_nx{
#
#====================================================================================
#
# This subroutine prints the Newton-X header followed by the content of 
# '$NX/version.text' file.
#
#------------------------------------------------------------------------------------
#
    ($mld)=@_;
#
    print_STDOUT("\n");
    print_STDOUT("          ============================================================  \n");
    print_STDOUT("                                     NEWTON-X                           \n");
    print_STDOUT("                   Newtonian dynamics close to the crossing seam        \n");
    print_STDOUT("          ============================================================  \n\n");
    print_STDOUT("                                  www.newtonx.org                       \n\n");
#
    open(VS, "$mld/version.txt") or die "Cannor open $mld/version.txt!";
    while(<VS>){print_STDOUT($_);}
    close(VS);
    
    print_STDOUT("\nNewton-X path: $mld\n\n");
}

sub end_nx{
#
#====================================================================================
#
# This subroutine end Newton-X. 
#
#------------------------------------------------------------------------------------
#
  print_STDOUT(
"          ====================== NEWTON-X ends here ==================  \n");
}

sub print_STDOUT{
#
#====================================================================================
#
# This subroutine prints the content of input variable $text in the Newton-X output 
# files, e.g., initcond.log, moldyn.log, ../RESULTS/nx.log, etc... according to 
# whether $text is the only input variable and/or the writing flag variable is on. 
#
#------------------------------------------------------------------------------------
#
    my ($istep, $kt, $flag, $eps, $text, $log_file, $n);
    $eps = 1E-9;
    ($text, $istep, $kt) = @_;
    if ( -f "auxp"){
       # .auxp is written by moldyn.pl flag=0 means print flag>0 means no print
       open(AP, "auxp") or die "cannot read auxp, ERROR in print_STDOUT\n$!\n";
       $_ = <AP>;
       chomp;
       s/\s*//g;
       $flag = $_;
    }else{
       # .if there is no auxp print anyway
       $flag = 0;
    }
    #. Variable $n stores the number of input variables actually introduced
    $n = @_;
    if (($n == 1) or ($n == 3)){
       # ..if only text is given, print always
       if ($n == 1){
          $flag = 0;
       }
       # ..$flag = $istep % $kt;
       if ($flag < $eps){
          print STDOUT "$text";
          if (-e "../RESULTS"){
             $log_file = "../RESULTS/nx.log";
             open(SECOUT,">>$log_file");
             print SECOUT "$text";
             close(SECOUT);
          }
       }
    }else{
       die "ERROR: wrong number of variables in print_STDOUT.\n";
    }
}

sub printf_STDOUT{
#
#====================================================================================
#
# Print formated outputs to moldyn.log and RESULTS/nx.log.
# Usage:
# printf_STDOUT($istep,$kt,$format,@output); print only if timestep ISTEP is multiple of KT
#
#------------------------------------------------------------------------------------
#
  my ($istep,$kt,$flag,$eps,$text,$log_file);
  $eps=1E-9;
  ($istep,$kt,$format,@text)=@_;
  if ( -f "auxp"){
    # auxp is written by moldyn.pl flag=0 means print flag>0 means no print
    open(AP,"auxp") or die "cannot read auxp, ERROR in print_STDOUT\n$!\n";
    $_=<AP>; chomp; s/\s*//g; $flag=$_;
  }else{
    # if there is no auxp print anyway
    $flag=0;
  }
  #  $flag=$istep % $kt;
  if ($flag < $eps){
    printf STDOUT $format,@text;
    if (-e "../RESULTS"){
      $log_file="../RESULTS/nx.log";
      open(SECOUT,">>$log_file");
      printf SECOUT $format,@text;
      close(SECOUT);
    }
  }
}

sub exit_error{
#
#====================================================================================
#
# This subroutine exits the program with an error message ($message). 
#
#------------------------------------------------------------------------------------
#
    local($mdle, $message,);
    $message = $_[0];
    $mdle    = $_[1];
    print_STDOUT("\n*******************************************************************************\n");
    print_STDOUT("     $mdle <ERROR> $message \n");
    print_STDOUT("*******************************************************************************\n");
    die "\n ERROR: NX is dying now ... \n";
}

sub err_message{
#
#====================================================================================
#
# Write error messages to runnx.error.
#
#------------------------------------------------------------------------------------
#
   local($scrv,$cbus);
   $cbus= $ENV{"COLUMBUS"};
   $scrv=0;
   $inderr = "inderr";
   if (-s $inderr){
      open(IE,"inderr") or die "Cannot open $inderr file!";
      $_=<IE>;
      chomp;
      $scrv = $_;
   }
   open(NE,">>$BASEDIR/../DEBUG/runnx.error");
   if ($scrv == 666){    # in moldyn03
      print NE
"Total energy changed more than Etot_jump in 1 timestep.
Options to solve (or avoid) this kind of problem are:
1 - reduce timestep dt;
2 - increase Etot_jump;
3 - check whether the electronic structure method and level
are adequate to compute the energies for the current geometry. \n";
   }
    if ($scrv == 667){    # in moldyn03
      print NE
"Total energy changed more than Etot_drift since t = 0.
Options to solve (or avoid) this kind of problem are:
1 - reduce timestep dt;
2 - increase Etot_drift;
3 - check whether the electronic structure method and level
are adequate to compute the energies for the current geometry. \n";
   }
    if ($scrv == 668){    # in sh
      print NE
"Total adiabatic population deviates more than popdev from the unity.
This means something is wrong with the integration of the TDSE.
Options to solve (or avoid) this kind of problem are:
1 - reduce timestep dt;
2 - increase ms;
3 - increase popdev;
4 - check whether the electronic structure method and level
are adequate to compute the energies, gradients and nonadiabatic
coupling vectors for the current geometry. \n";
   }
    if ($scrv == 669){    # in sh_mmod
      print NE
"Hopping probability is larger than one (see file RESULTS/tprob).
This means something is wrong with the integration of the TDSE or
with the nonadiabatic coupling vectors. (Is the system in a conical
intersectio? Surface hopping may fail there.)
Options to solve (or avoid) this kind of problem are:
1 - reduce timestep dt;
2 - increase ms;
3 - try another TDSE integration method (change integrator);
4 - check whether the electronic structure method and level
are adequate to compute the energies, gradients and nonadiabatic
coupling vectors for the current geometry. \n";
   }
    if ($scrv == 670){    # in sh_mmod
      print NE
"Hopping probability is lower than zero (see file RESULTS/tprob).
This means something is wrong with the integration of the TDSE or
with the nonadiabatic coupling vectors. (Is the system in a conical
intersection? Surface hopping may fail there.)
Options to solve (or avoid) this kind of problem are:
1 - reduce timestep dt;
2 - increase ms;
3 - try another TDSE integration method (change integrator);
4 - check whether the electronic structure method and level
are adequate to compute the energies, gradients and nonadiabatic
coupling vectors for the current geometry. \n";
   }
    if ($scrv == 671){
      print NE
"The COLUMBUS version ($cbus)
might not be compatible with the current version of Newton-X
for this specific kind of job (ciudg in control.run and
NSTAT > NSTATDYN in control.dyn).
If in the first step of this job COLUMBUS ends with error message:
\"error in cigrd: eci and eci(eff) do not match\"
You should contact the Columbus distributors (www.univie.ac.at/columbus)
and require an updated version of the runc file.
\n";
   }
   close(NE);
   close(IE);
}

sub open_err{
#
#====================================================================================
#
# This subroutine is called when opening a file leads to an error. Use it as:
# 
# open_err(filename, intent, warn/die, module/routine)
#
# - place "warn" as third argument if you want the program to continue.
# - if you want to omit an argument, give "" (e.g. open_err("","","","")) 
# - you can omit argument from the end completely (open_err("filename") is valid). 
# Standard values will be inserted: 
#         filename = "file"
#         intent = "unkown intent" 
#         warn/die = "die" 
#         mod/routine = "" 
#
#------------------------------------------------------------------------------------
#
     my ($filename, $intent, $action, $cmdle); #read/write,warn/die,calling mod/routine
#
     if ($#_ == 3){
        ($filename, $intent, $action, $cmdle) = @_;
     }elsif ($#_ == 2){
        $cmdle = "";
        ($filename, $intent, $action) = @_;
     }elsif ($#_ == 1){
        $cmdle = "";
        $action = "";
        ($filename, $intent) = @_;
     }elsif ($#_ == 0){
        $cmdle = "";
        $action = "";
        $intent = "";
       ($filename) = @_;
     }
#
     if ($filename eq ""){$filename = "file";}
     if ($intent eq ""){$intent = "unknown intent";}
#
     if ($action eq "warn"){
        print_STDOUT("$cmdle WARNING: could not open $filename for $intent\n");
     }else{
        die "$cmdle could not open $filename for $intent\n";
     }
 }

sub close_err{
#
#====================================================================================
#
# This subroutine is called when closing a file leads to an error.
# Usage:  close_err(filename)
#
#------------------------------------------------------------------------------------
#
    local ($filename);
    die "Could not close $filename\nSystems answer was: $!";
}

#
#====================================================================================
#       FILES AND DIRECTORIES MANIPULATION
#====================================================================================
#

sub rm_files{
#
#====================================================================================
#
# Delete files
# usage: rm_files(path,type,size)
# path - directory path
# type - B for binary (case insensitive)
#        T for text
#        A for all
# size - remove only files above this size in MB
# Example:
# rm_files("TEMP/WORK","B",1)
# will remove all binary files from TEMP/WORK whose size is
# larger or equal to 1 MB.
# Mario Barbatti, Mar 2007
#
#------------------------------------------------------------------------------------
#
   local (@files);
   local ($element,$dir_path,$fname,$max_size,$MB2byte,$type,$value);
   ($dir_path,$type,$value)=@_;
   $type = lc $type;

   if (!-e $dir_path){
     return;
   }
   $MB2byte = 1048576;
   $max_size = $value*$MB2byte;

   opendir (DIR,$dir_path) or die "Can't open current dir: $dir_path\n";
   @files = grep (!/^\.\.?$/, readdir (DIR));
   closedir (DIR);

   for $element (@files) {
       $fname = "$dir_path/$element";
       if (-B $fname){
          if (-s $fname >= $max_size){
             if (($type eq "b") or ($type eq "a")){
               system("rm -rf $fname");   # Delete binary files
             }
          }
       }elsif (-T $fname){
          if (-s $fname >= $max_size){
             if (($type eq "t") or ($type eq "a")){
               system("rm -rf $fname");   # Delete text files
             }
          }
       }
   }
}

sub save_files{
#
#====================================================================================
#
# This subroutine saves to SAVE_TEMP directory a list of files defined in 
# 'save_file' file.
#
#------------------------------------------------------------------------------------
#
     my ($file, $save, @g, $last);
     $save = "save_file";
     $ST = "SAVE_TEMP";
     if (-s $save){
        if (!-s $ST){
           system("mkdir $ST");
        }
        $i = 0;
        open(SF, $save) or die "Cannot read $save";
        while(<SF>){
             chomp;
             $_ =~ s/^\s*//;
             $_ =~ s/\s*$//;
             $file = $_;
             if (-s $file){
                system("cp -rf $file $ST/.");
                if ($file =~ /\//){
                   @g = split(/\//, $file);
                   $last = $#g;
                   $file = $g[$last];
                }
                system("cd $ST; tar --remove-files -zcf $file.tgz $file");
             }
        }
        close(SF);
     }
 }


sub cptemp{
#
#====================================================================================
#
# This subroutine tries to save TEMP before dying, but only once.
#
#------------------------------------------------------------------------------------
#
  my (@status);
  $ctd = "control.d";
  if (!-e "savetemp") {
    $t = "";
    if (-e $ctd) {
      @status=load_status($ctd,$mdle);
      $istep = $status[1];
      $t     = $status[7];
      $t     = $t*1;
      print_STDOUT("\n$mdle   ::ERROR::   step = $istep time = $t fs \n");
    }
    if (-e "../DEBUG/COL.$t") {
      $direct = "../DEBUG/COL.$t";
      rm_files("WORK","t",5); # Delete text files above 5 MB
      rm_files("WORK","b",0); # Delete all binary files
    } elsif(-e "../DEBUG/G03.$t") {
      $direct = "../DEBUG/G03.$t";
    } elsif(-e "../DEBUG/G09.$t") {
      $direct = "../DEBUG/G09.$t";
    } else {
      $direct = "../DEBUG";
      rm_files("../TEMP","a",5); # Delete all binary files
    }
    print_STDOUT("Trying to save TEMP directory to $direct \n");
    system("cp -rf ../TEMP $direct/.");
    system("touch savetemp");
  }
}

sub saveoldfile{
#
#====================================================================================
#
# This subroutine copies the contents of an input file to a temp file 'tmp-clp'. 
# If the old file is not present an empty tmp-clp is created. 
#
#------------------------------------------------------------------------------------
#
     local ($file);
     $file = $_[0];
     $/ = "\n";
     if (open OLD, "<$file"){
        #print_STDOUT("file $file present\n");
        open TMP, ">tmp-clp" or die "cannot create temporary file\n sys answer was: $!";
        while(<OLD>){
             #print_STDOUT($_);
             print TMP $_;
        }
        close OLD or die "could not close $file\nsys answer was: $!";
        close TMP or die "could not close temporary file\n sys answer was: $!";
     }else{
        #print_STDOUT("no file $file - creating a new keyfile\n)";
        system ("rm -f tmp-clp");
        system ("touch tmp-clp");
     }
 }

#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#====================================================================================
#
#                    THIRD-PARTY PROGRAMS (TPP) SECTION
#
#         The subroutines defined below are related to the TPP programs.
#
#====================================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#

#
#====================================================================================
#       GENERAL TPP ROUTINES
#====================================================================================
#

sub prog_config{
#
#====================================================================================
#
# This subroutine returns third-party program configuration according to $prog. 
# Use it as: %progconf = prog_config($prog);
#
# For each program:
# $progconf{progname}  = string to identify program/ during NX execution.
# $progconf{methodname}= string to identify method during NX execution.
# $progconf{ic}        = y/n: display in initial conditions menu of nxinp?
# $progconf{dyn}       = y/n: display in dynamics menu of nxinp?
# $progconf{hyb}       = y/n: display in hybrid calculations menu of nxinp?
# $progconf{ip}        = y/n: display in photoelectron spectrum menu of nxinp?
# $progconf{key}       = numerical value associated to the program.
# $progconf{method}    = string specifying method.
# $progconf{parfile}   = name of the parameter file .par.
# $progconf{nad_exec}  = y/n: enabled for nonadiabatic dynamics (TD-BA included).
# $progconf{lvprt_d}   = default value of lvprt keyword.
# $progconf{vdoth_d}   = default value of vdoth keyword.
# $progconf{never_state_d} = default value of never_state keyword.
# $progconf{cio_options_d} = default value of cio_options keyword. 
# $progconf{cisc_options_d} = default value of cisc_options keyword.
# $progconf{cprog_d}   = default value of cprog keyword (cioverlap program).
# $progconf{progic}    = program for computing initial conditions.
# $progconf{progdyn}   = program for computing dynamics.
#
#------------------------------------------------------------------------------------
#
    my ($prog,%progconf,$progmax);

    ($prog)  = @_;

    if (($prog >=-0.05) and ($prog < 0.05)){
      %progconf=(progname   => "analytical",
	         methodname => "NULL",
	         ic         => "n",
                 dyn        => "y",
                 hyb        => "y",
                 ip         => "n",
                 key        => sprintf("%4d",0),
                 label      => "ANALYTICAL",
		 method     => "NULL",
		 parfile    => "analyt.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 0,
                 never_state_d  => 0,
                 cio_options_d  => "NULL", 
                 cisc_options_d => "NULL",
                 cprog_d    => 0,
	         progic     => "NULL",
	         progdyn    => "run-analyt.pl");
    }
    if (($prog >= 0.95) and ($prog < 1.05)){
      %progconf=(progname   => "columbus",
	         methodname => "columbus_sa-mcscf",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "y",
                 ip         => "n",
                 key        => sprintf("%4.1f",1.0),
                 label      => "COLUMBUS",
		 method     => "SA-MCSCF GRADIENT",
		 parfile    => "columbus.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 0,
                 never_state_d  => 0,
                 cio_options_d  => "\"-b -t 5e-4 -e 2 -i\"", 
                 cisc_options_d => "\"\"",
                 cprog_d    => 1,
                 progic     => "run_col_initcond.pl",
                 progdyn    => "run-col.pl");	    
    }
    if (($prog >= 1.05) and ($prog < 1.15)){
      %progconf=(progname   => "columbus",
                 methodname => "columbus_mrci",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "y",
                 ip         => "n",
                 key        => sprintf("%4.1f",1.1),
                 label      => "COLUMBUS",
                 method     => "MRCI GRADIENT",
                 parfile    => "columbus.par",
                 nad_exec   => "y",
		 lvprt_d    => 1,
                 vdoth_d    => 0,
                 never_state_d  => 0,
                 cio_options_d  => "\"-b -t 5e-4 -e 2 -i\"",
                 cisc_options_d => "\"\"",
                 cprog_d    => 1,
                 progic     => "run_col_initcond.pl",
                 progdyn    => "run-col.pl");
    }
    if (($prog >= 1.95) and ($prog < 2.05)){
      %progconf=(progname   => "turbomole",
	         methodname => "turbomole-ricc2",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "y",
                 ip         => "n",
                 key        => sprintf("%4.1f",2.0),
                 label      => "TURBOMOLE",
		 method     => "CC2",
		 parfile    => "turbomole.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 1,
                 never_state_d  => 1,
                 cio_options_d  => "\"-s transmomin -a -t 5e-4 -e -1\"", 
                 cisc_options_d => "\"\"",
                 cprog_d    => 2,
                 progic     => "run_tm_initcond.pl",
                 progdyn    => "run-turbo.pl");
    }
    if (($prog >= 2.05) and ($prog < 2.15)){
      %progconf=(progname   => "turbomole",
	         methodname => "turbomole-tddft",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "y",
                 ip         => "n",
                 key        => sprintf("%4.1f",2.1),
                 label      => "TURBOMOLE",
		 method     => "TDDFT",
		 parfile    => "turbomole.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 1,
                 never_state_d  => 1,
                 cio_options_d  => "\"-s transmomin -a -t 5e-4 -e -1\"", 
                 cisc_options_d => "\"-o\"",
                 cprog_d    => 2,
                 progic     => "run_tm_initcond.pl",
                 progdyn    => "run-turbo.pl");
    }
    if (($prog >= 2.15) and ($prog < 2.25)){
      %progconf=(progname   => "turbomole",
	         methodname => "turbomole-riadc2",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "y",
                 ip         => "n",
                 key        => sprintf("%4.1f",2.2),
                 label      => "TURBOMOLE",
		 method     => "ADC(2)",
		 parfile    => "turbomole.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 1,
                 never_state_d  => 1,
                 cio_options_d  => "\"-s transmomin -a -t 5e-4 -e -1\"", 
                 cisc_options_d => "\"\"",
                 cprog_d    => 2,
                 progic     => "run_tm_initcond.pl",
                 progdyn    => "run-turbo.pl");
    }
    if (($prog >= 2.95) and ($prog < 3.05)){
      %progconf=(progname   => "aces2",
	         methodname => "NULL",
                 ic         => "n",
                 dyn        => "n",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4d",3),
                 label      => "ACES II",
		 method     => "NULL",
		 parfile    => "aces2.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 0,
                 never_state_d  => 0,
                 cio_options_d  => "NULL", 
                 cisc_options_d => "NULL",
                 cprog_d    => 0,
                 progic     => "NULL",
                 progdyn    => "run-a2.pl");
    }
    if (($prog >= 3.95) and ($prog < 4.05)){
      %progconf=(progname   => "mopac",
	         methodname => "FOMO-CI",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4d",4),
                 label      => "MOPAC",
		 method     => "FOMO-CI",
		 parfile    => "mopac.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => -1,
                 never_state_d  => 0,
                 cio_options_d  => "NULL", 
                 cisc_options_d => "NULL",
                 cprog_d    => 3,
                 progic     => "run_mopac_initcond.pl",
                 progdyn    => "run-mopac.pl");
    }
    if (($prog >= 4.45) and ($prog < 4.55)){
      %progconf=(progname   => "exc_mopac",
	         methodname => "EXASH FOMO-CI/TINKER",
                 ic         => "n",
                 dyn        => "y",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4.1f",4.5),
                 label      => "MOPAC/TINKER",
		 method     => "EXASH",
		 parfile    => "exc.inp",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => -1,
                 never_state_d  => 0,
                 cio_options_d  => "NULL", 
                 cisc_options_d => "NULL",
                 cprog_d    => 3,
                 progic     => "NULL",
                 progdyn    => "run-exc_mopac.pl");
    }
    if (($prog >= 4.95) and ($prog < 5.05)){
      %progconf=(progname   => "orca",
	         methodname => "orca-tddft",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4d",5),
                 label      => "ORCA",
		 method     => "TDDFT",
		 parfile    => "orca.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 2,
                 never_state_d  => 0,
                 cio_options_d  => "NULL", 
                 cisc_options_d => "NULL",
                 cprog_d    => 0,
                 progic     => "run_orca_initcond.pl",
                 progdyn    => "run-orca.pl");
    }
    if (($prog >= 5.95) and ($prog < 6.05)){
      %progconf=(progname   => "gau",
	         methodname => "gaussian-cas",
                 ic         => "n",
                 dyn        => "y",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4.1f",6.0),
                 label      => "GAUSSIAN",
		 method     => "CASSCF",
		 parfile    => "gau.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 0,
                 never_state_d  => 0,
                 cio_options_d  => "NULL", 
                 cisc_options_d => "NULL",
                 cprog_d    => 0,
                 progic     => "NULL",
                 progdyn    => "run-gau-cas.pl");
    }
    if (($prog >= 6.45) and ($prog < 6.55)){
      %progconf=(progname   => "gau",
	         methodname => "gaussian-lr",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "n",
                 ip         => "y",
                 key        => sprintf("%4.1f",6.5),
                 label      => "GAUSSIAN",
		 method     => "TDDFT",
		 parfile    => "gau.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 1,
                 never_state_d  => 1,
                 cio_options_d  => "\"-s transmomin -a -t 5e-4 -e -1\"", 
                 cisc_options_d => "\"-o\"",
                 cprog_d    => 2,
                 progic     => "run_gau_initcond.pl",
                 progdyn    => "run-gau.pl");
    }
    if (($prog >= 6.95) and ($prog < 7.05)){
      %progconf=(progname   => "tinker",
	         methodname => "NULL",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "y",
                 ip         => "n",
                 key        => sprintf("%4d",7),
                 label      => "TINKER",
                 method     => "MM",
		 parfile    => "tinker.par",
		 nad_exec   => "n",
		 lvprt_d    => 1,
		 vdoth_d    => 2,
                 never_state_d  => 0,
                 cio_options_d  => "NULL", 
                 cisc_options_d => "NULL",
                 cprog_d    => 0,
                 progic     => "run_tinker_initcond.pl",
                 progdyn    => "run-tinker.pl");
    }
    if (($prog >= 7.95) and ($prog < 8.05)){
      %progconf=(progname   => "NULL",
	         methodname => "NULL",
                 ic         => "n",
                 dyn        => "n",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4.1f",8.0),
                 label      => "NULL",
		 method     => "NULL",
		 parfile    => "NULL",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 2,
                 never_state_d  => 0,
                 cio_options_d  => "NULL", 
                 cisc_options_d => "NULL",
                 cprog_d    => 0,
                 progic     => "NULL",
                 progdyn    => "NULL");
    }
    if (($prog >= 8.05) and ($prog < 8.15)){
      %progconf=(progname   => "NULL",
	         methodname => "NULL",
                 ic         => "n",
                 dyn        => "n",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4.1f",8.1),
                 label      => "NULL",
		 method     => "NULL",
		 parfile    => "NULL",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 2,
                 never_state_d  => 0,
                 cio_options_d  => "NULL", 
                 cisc_options_d => "NULL",
                 cprog_d    => 0,
                 progic     => "NULL",
                 progdyn    => "NULL");
    }
    if (($prog >= 8.45) and ($prog < 8.55)){
      %progconf=(progname   => "dftb+",
	         methodname => "NULL",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4.1f",8.5),
                 label      => "DFTB+",
		 method     => "TD-DFTB",
		 parfile    => "dftb+.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 1,
                 never_state_d  => 1,
                 cio_options_d  => "\"-s transmomin -a -t 5e-3 -e -1\"", 
                 cisc_options_d => "\"-o\"",
                 cprog_d    => 2,
                 progic     => "run_dftb+_initcond.pl",
                 progdyn    => "run-dftb+.pl");
    }
    if (($prog >= 8.95) and ($prog < 9.05)){
      %progconf=(progname   => "dftci",
	         methodname => "NULL",
                 ic         => "y",
                 dyn        => "n",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4d",9),
                 label      => "DFTCI",
		 method     => "DFT/MRCI",
		 parfile    => "dftci.par",
		 nad_exec   => "n",
		 lvprt_d    => 1,
		 vdoth_d    => 2,
                 never_state_d  => 0,
                 cio_options_d  => "NULL", 
                 cisc_options_d => "NULL",
                 cprog_d    => 0,
                 progic     => "run_dftci_initcond.pl",
                 progdyn    => "NULL");
    }
    if (($prog >= 9.95) and ($prog < 10.05)){
      %progconf=(progname   => "gamess",
	         methodname => "GAMESS-MCSCF-NAD",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4.1f",10.0),
                 label      => "GAMESS",
		 method     => "MCSCF (NONADIABATIC)",
		 parfile    => "gamess.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 0,
                 never_state_d  => 0,
                 cio_options_d  => "NULL", 
                 cisc_options_d => "NULL",
                 cprog_d    => 0,
                 progic     => "run_gamess_initcond.pl",
                 progdyn    => "run-gamess.pl");
    }
    if (($prog >= 10.05) and ($prog < 10.15)){
      %progconf=(progname   => "gamess",
	         methodname => "GAMESS-ARB-AD",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4.1f",10.1),
                 label      => "GAMESS",
		 method     => "ARBITARY METHOD (ADIABATIC)",
		 parfile    => "gamess.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 2,
                 never_state_d  => 0,
                 cio_options_d  => "NULL", 
                 cisc_options_d => "NULL",
                 cprog_d    => 0,
                 progic     => "run_gamess_initcond.pl",
                 progdyn    => "run-gamess.pl");
    }
    if (($prog >= 10.95) and ($prog < 11.05)){
      %progconf=(progname   => "mlatom",
	         methodname => "NULL",
                 ic         => "n",
                 dyn        => "y",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4d",11),
                 label      => "MLATOM",
                 method     => "ML MODEL",
		 parfile    => "mlatom.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 0,
                 never_state_d  => 0,
                 cio_options_d  => "NULL", 
                 cisc_options_d => "NULL",
                 cprog_d    => 0,
                 progic     => "NULL",
                 progdyn    => "run-mlatom.pl");
    }
    if (($prog >= 11.95) and ($prog < 12.05)){
      %progconf=(progname   => "bagel",
	         methodname => "BAGEL-MR-NAD",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "y",
                 ip         => "n",
                 key        => sprintf("%4.1f",12.0),
                 label      => "BAGEL",
		 method     => "MCSCF/CASPT2 (NONADIABATIC)",
		 parfile    => "bagel.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 0,
                 never_state_d  => 0,
                 cio_options_d  => "NULL", 
                 cisc_options_d => "NULL",
                 cprog_d    => 0,
                 progic     => "run_bagel_initcond.pl",
                 progdyn    => "run-bagel.pl");
    }
    if (($prog >= 12.05) and ($prog < 12.15)){
      %progconf=(progname   => "bagel",
	         methodname => "BAGEL-ARB-AD",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "y",
                 ip         => "n",
                 key        => sprintf("%4.1f",12.1),
                 label      => "BAGEL",
		 method     => "ARBITRARY METHOD (ADIABATIC)",
		 parfile    => "bagel.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 2,
                 never_state_d  => 0,
                 cio_options_d  => "NULL", 
                 cisc_options_d => "NULL",
                 cprog_d    => 0,
                 progic     => "run_bagel_initcond.pl",
                 progdyn    => "run-bagel.pl");
    }
    if (($prog >= 12.95) and ($prog < 13.05)){
      %progconf=(progname   => "pysoc",
                 methodname => "NULL",
                 ic         => "y",
                 dyn        => "n",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4d",13),
                 label      => "PYSOC",
                 method     => "NULL",
                 parfile    => "pysoc.par",
                 nad_exec   => "y",
		 lvprt_d    => 2,
                 vdoth_d    => 2,
                 never_state_d  => 0,
                 cio_options_d  => "\"\"",
                 cisc_options_d => "\"\"",
                 cprog_d    => 1,
                 progic     => "run_pysoc_initcond.py",
                 progdyn    => "NULL");
    }
    if (($prog >= 13.95) and ($prog < 14.05)){
      %progconf=(progname   => "cp2k",
                 methodname => "NULL",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4d",14),
                 label      => "CP2K",
                 method     => "TDDFT",
                 parfile    => "cp2k.par",
                 nad_exec   => "y",
                 lvprt_d    => 1,
                 vdoth_d    => 1,
                 never_state_d  => 1,
                 cio_options_d  => "\"-s transmomin -a -t 5e-4 -e -1\"",
                 cisc_options_d => "\"-o -c\"",
                 cprog_d    => 2,
                 progic     => "run_cp2k_initcond.pl",
                 progdyn    => "run-cp2k.pl");
    }
    if (($prog >= 14.95) and ($prog < 15.05)){
      %progconf=(progname   => "turbomole",
                 methodname => "turbomole-ricc2",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4.1f",15.0),
                 label      => "FROMAGE(turbomole)",
                 method     => "CC2",
                 parfile    => "turbomole.par",
                 nad_exec   => "y",
                 lvprt_d    => 1,
                 vdoth_d    => 1,
                 never_state_d  => 1,
                 cio_options_d  => "\"-s transmomin -a -t 5e-4 -e -1\"",
                 cisc_options_d => "\"\"",
                 cprog_d    => 2,
                 progic     => "run_fromage_initcond.pl",
                 progdyn    => "run-fromage.pl");
    }
    if (($prog >= 15.05) and ($prog < 15.15)){
      %progconf=(progname   => "turbomole",
                 methodname => "turbomole-tddft",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4.1f",15.1),
                 label      => "FROMAGE(turbomole)",
                 method     => "TDDFT",
                 parfile    => "turbomole.par",
                 nad_exec   => "y",
                 lvprt_d    => 1,
                 vdoth_d    => 1,
                 never_state_d  => 0,
                 cio_options_d  => "\"-s transmomin -a -t 5e-4 -e -1\"",
                 cisc_options_d => "\"-o\"",
                 cprog_d    => 2,
                 progic     => "run_fromage_initcond.pl",
                 progdyn    => "run-fromage.pl");
    }
    if (($prog >= 15.15) and ($prog < 15.25)){
      %progconf=(progname   => "turbomole",
                 methodname => "turbomole-riadc2",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4.1f",15.2),
                 label      => "FROMAGE(turbomole)",
                 method     => "ADC(2)",
                 parfile    => "turbomole.par",
                 nad_exec   => "y",
                 lvprt_d    => 1,
                 vdoth_d    => 1,
                 never_state_d  => 1,
                 cio_options_d  => "\"-s transmomin -a -t 5e-4 -e -1\"",
                 cisc_options_d => "\"\"",
                 cprog_d    => 2,
                 progic     => "run_fromage_initcond.pl",
                 progdyn    => "run-fromage.pl");
    }
    if (($prog >= 15.35) and ($prog < 15.45)){
      %progconf=(progname   => "gau",
                 methodname => "gaussian-lr",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4.1f",15.4),
                 label      => "FROMAGE(G16)",
                 method     => "TDDFT",
                 parfile    => "gau.par",
                 nad_exec   => "y",
                 lvprt_d    => 1,
                 vdoth_d    => 1,
                 never_state_d  => 1,
                 cio_options_d  => "\"-s transmomin -a -t 5e-4 -e -1\"",
                 cisc_options_d => "\"-o\"",
                 cprog_d    => 2,
                 progic     => "run_fromage_initcond.pl",
                 progdyn    => "run-fromage.pl");
    }
    if (($prog >= 15.45) and ($prog < 15.55)){
      %progconf=(progname   => "openmolcas",
                 methodname => "openmolcas-cas",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4.1f",15.5),
                 label      => "FROMAGE(OMolcas)",
                 method     => "CASSCF-CASPT2",
                 parfile    => "openmolcas.par",
                 nad_exec   => "y",
                 lvprt_d    => 1,
                 vdoth_d    => 0,
                 never_state_d  => 0,
                 cio_options_d  => "NULL",
                 cisc_options_d => "NULL",
                 cprog_d    => 3,
                 progic     => "run_fromage_initcond.pl",
                 progdyn    => "run-fromage.pl");
    }
    if (($prog >= 15.55) and ($prog < 15.65)){
      %progconf=(progname   => "orca",
                 methodname => "orca-tddft",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "n",
                 ip         => "n",
                 key        => sprintf("%4.1f",15.6),
                 label      => "FROMAGE(Orca)",
                 method     => "TDDFT",
                 parfile    => "orca.par",
                 nad_exec   => "y",
                 lvprt_d    => 1,
                 vdoth_d    => 2,
                 never_state_d  => 0,
                 cio_options_d  => "NULL",
                 cisc_options_d => "NULL",
                 cprog_d    => 3,
                 progic     => "run_fromage_initcond.pl",
                 progdyn    => "run-fromage.pl");
    }
    if (($prog >= 19.95) and ($prog < 20.05)){
      %progconf=(progname   => "hybrid",
	         methodname => "NULL",
                 ic         => "y",
                 dyn        => "y",
                 hyb        => "y",
                 ip         => "n",
                 key        => sprintf("%4d",20),
                 label      => "HYBRID",
		 method     => "NULL",
		 parfile    => "hybrid.par",
		 nad_exec   => "y",
		 lvprt_d    => 1,
		 vdoth_d    => 2,
                 never_state_d  => 0,
                 cio_options_d  => "\"\"", 
                 cisc_options_d => "\"\"",
                 cprog_d    => 1,
                 progic     => "run_hybrid_initcond.pl",
                 progdyn    => "run-hybrid.pl");
    }
    return %progconf;
}

sub check_thirdparty{
#
#====================================================================================
# This subroutine checks the third-party inputs
# For each external program, we have a block like:
# if ($progname eq <ExtProg>){
#   if (ExtProg has defined parameters){ 
#     . Load default parameter values for ExtProg
#     . Display parameter values in NX output
#     . Rewrite parameter file 
#   }   
#   Additionally: 
#   . Check ExtProg environment variables [optional]
#   . Check input consistency [optional] 
#   . Print credits for theinterface [optional]   
# }
#
# This subroutine is called by: 
#   . moldyn.pl ($type = "dyn")
#   . initcond.pl ($type = "ic")
#   . nxinp ($type = "mic" for initial conditions menu and $type = "mdy" for dynamics menu).
#
     my(@line);
     my ($path,$am1_file,$am1_alpha,$am1_beta,$am1_kx,$am1_ky);
     my ($am1_delta,$am1_x1,$am1_x2,$am1_x3,$am1_gamma);
     my (@g,$g_vers);
     my ($cbus,$gau,$gaur,$mlatom,$dftbp,$value_found,$value_test,@control);
     my ($type,$job,$listpar,$progpar,$prog);
     my ($BASEDIR,$TP,$JAD,$nstat,$mdle,$lvprt,$thres);
     my (%progconf,$hybrid_tpp);
     my ($cp2k);

     $TP="TEMP";
     $JAD="JOB_AD";
     $JND="JOB_NAD";

     ($type,$prog,$BASEDIR,$hybrid_tpp,$nstat,$mdle,$lvprt,$thres,$anmod,$tully_mod,$cs_mod,$vdoth)=@_;
     %hybrid_tpp = %{$hybrid_tpp};

     if ($prog eq ""){exit_error("PROG = $prog is out of range and does not correspond to any program.",$mdle);}

     if ($type eq "dyn"){
       $job = "Dynamics";
     }
     if ($type eq "ic"){
       $job = "Initial conditions";
     }

     %progconf   = prog_config($prog);
     $progname   = $progconf{progname};
     $methodname = $progconf{methodname};
     $label      = $progconf{label};
     $progpar    = $progconf{parfile};

     if (($type eq "dyn") or ($type eq "ic")){
       print_STDOUT("$job with $label\n\n");
       $credits = read_credits("$progname");
       if (defined $credits){
         print_STDOUT("$credits");
       }
     }

     ## ANALYTICAL ##
     #----------------------------------------------------------------------------
     if (($progname eq "analytical") or ($hybrid_tpp{analytical} == 1)) {
       
       if ($hybrid_tpp{analytical} == 1){
	  $prog = 0;
          %progconf   = prog_config($prog);
          $progpar    = $progconf{parfile};
          if (-s "$BASEDIR/$progpar"){system("cp -f $BASEDIR/$progpar .");}
       }

       # Load parameters
       ($anmod_d,$path_d,$tully_mod_d,$cs_mod_d)=load_defaults_tpp($prog,"main","");
       $path      = getkeyword($progpar,"path"      ,$path_d);
       if (($type eq "ic") or ($type eq "dyn")){
         $anmod     = getkeyword($progpar,"anmod"     ,$anmod_d);
         $tully_mod = getkeyword($progpar,"tully_mod" ,$tully_mod_d);
         $cs_mod    = getkeyword($progpar,"cs_mod"    ,$cs_mod_d);
       }

       # Parameter list
       $listpar = " anmod      = $anmod\n"
                 ." path       = $path\n";

       if ($anmod eq "tully_models.pl"){
         $listpar = $listpar
                   ." tully_mod  = $tully_mod\n";
       }elsif($anmod eq "model_cs_fssh"){
         $listpar = $listpar
                   ." cs_mod     = $cs_mod\n"; 
       }

       # Display
       if (($type eq "ic") or ($type eq "dyn")){
         print_STDOUT("$listpar\n");
       }

       # Write new parameter file
       open(INP,">$progpar") or die "$mdle Cannot write $progpar.";
       print INP "$listpar";
       close(INP);

     }

     ## COLUMBUS ##
     #----------------------------------------------------------------------------
     if (($progname eq "columbus") or ($hybrid_tpp{columbus} == 1)) {

       if ($hybrid_tpp{columbus} == 1){
	  $prog = 1;
          %progconf   = prog_config($prog);
          $progpar    = $progconf{parfile};
          if (-s "$BASEDIR/$progpar"){system("cp -f $BASEDIR/$progpar .");}
       }

       # Check environment variable
       if (($type eq "dyn") or ($type eq "ic")){
         $cbus  = $ENV{"COLUMBUS"};  # Columbus environment
         if (!-e "$cbus/runc"){exit_error("Cannot find Columbus program. Check \$COLUMBUS variable.",$mdle);}
       }

       # Load parameters
       ($memcol_d,$cirestart_d,$reduce_tol_d,$quad_conv_d,$mc_conv_d,$ci_conv_d,$ivmode_d,$mocoef_d,$prt_mo_d,$all_grads_d,$grad_lvl_d,$cver_d,$ci_iter_d,$citol_d,$mull_pop_d)=load_defaults_tpp($prog,"","");
       $memcol     = getkeyword($progpar,"mem"       ,$memcol_d);
       $cver       = getkeyword($progpar,"cver"      ,$cver_d);
       if (($type eq "dyn") or ($type eq "mdy")){
         $cirestart  = getkeyword($progpar,"cirestart" ,$cirestart_d);
         $reduce_tol = getkeyword($progpar,"reduce_tol",$reduce_tol_d);
         $quad_conv  = getkeyword($progpar,"quad_conv" ,$quad_conv_d);
         $mc_conv    = getkeyword($progpar,"mc_conv"   ,$mc_conv_d);
         $ci_conv    = getkeyword($progpar,"ci_conv"   ,$ci_conv_d);
         $ivmode     = getkeyword($progpar,"ivmode"    ,$ivmode_d);
         $mocoef     = getkeyword($progpar,"mocoef"    ,$mocoef_d);
         $prt_mo     = getkeyword($progpar,"prt_mo"    ,$prt_mo_d);
         $all_grads  = getkeyword($progpar,"all_grads" ,$all_grads_d);
	 $grad_lvl   = getkeyword($progpar,"grad_lvl"  ,$grad_lvl_d);
         $cver       = getkeyword($progpar,"cver"      ,$cver_d);
         $ci_iter    = getkeyword($progpar,"ci_iter"   ,$ci_iter_d);
         $citol      = getkeyword($progpar,"citol"     ,$citol_d);
         $mull_pop   = getkeyword($progpar,"mull_pop"  ,$mull_pop_d);
       }

       # Check input
       if ($type eq "dyn"){
         if ($all_grads == 1){
            if ($thres == 0){
               $all_grads = 0;
               print_STDOUT("$mdle For thres = 0, all_grads must be 0. Changing it.\n\n");
            }elsif($thres >0){
               if($lvprt >= 3){print_STDOUT("$mdle all_grads = 1: this option will not work for adiabatic jobs.\n\n");}
            }
         }
       }

       # Parameter list
       $listpar = " memcol     = $memcol\n"
                 ." cver       = $cver\n";

       if (($type eq "dyn") or ($type eq "mdy")){
         $listpar = $listpar
                   ." cirestart  = $cirestart\n"
                   ." quad_conv  = $quad_conv\n"
                   ." mc_conv    = $mc_conv\n"
                   ." ci_conv    = $ci_conv\n"
                   ." ivmode     = $ivmode\n"
                   ." mocoef     = $mocoef\n"
                   ." prt_mo     = $prt_mo\n"
                   ." ci_iter    = $ci_iter\n"
                   ." citol      = $citol\n"
                   ." reduce_tol = $reduce_tol\n"
                   ." all_grads  = $all_grads\n"
                   ." grad_lvl   = $grad_lvl\n";
       }

       # Display
       if (($type eq "ic") or ($type eq "dyn")){
         print_STDOUT("COLUMBUS path: $cbus\n");
         print_STDOUT("$listpar\n");
       }

       # Write new parameter file
       open(INP,">$progpar") or die "$mdle Cannot write $progpar.";
       print INP "$listpar";
       close(INP);

       # Check CIOVERLAP
       if (($vdoth == -1) or ($vdoth == 1)){
         my $mldcio = ciocheck($mld);
	 print_STDOUT("CIOVERLAP path: $mldcio");
       }

     }

     ## TURBOMOLE ##
     #----------------------------------------------------------------------------
     if (($progname eq "turbomole") or ($hybrid_tpp{turbomole} == 1)) {

       if ($hybrid_tpp{analytical} == 1){
	  $prog = 2;
          %progconf   = prog_config($prog);
          $progpar    = $progconf{parfile};
          if (-s "$BASEDIR/$progpar"){system("cp -f $BASEDIR/$progpar .");}
       }
      
       # Load parameters 
       ($parallel_d) = load_defaults_tpp($prog,"","");
       $parallel     = getkeyword($progpar,"parallel",$parallel_d); 

       # Parameter list
       $listpar = "parallel = $parallel\n";

       # Display
       if (($type eq "ic") or ($type eq "dyn")){
         print_STDOUT("$listpar\n");
       }

       # Write new parameter file
       open(INP,">$progpar") or die "$mdle Cannot write $progpar.";
       print INP "$listpar";
       close(INP);

       # Check input
       if (($type eq "ic") or ($type eq "dyn")){
         @control=();
         @control=read_control("$JAD");
         $value_found=get_control("ricc2",@control);
         print_STDOUT("Method: $methodname\n");
         $value_test="";
         if ($methodname eq "turbomole-ricc2"){
            $value_test="cc2";
         }elsif($methodname eq "turbomole-riadc2"){
            $value_test="adc";
         }
         if ($nstat > 1){ 
           if (($methodname eq "turbomole-ricc2") or ($methodname eq "turbomole-riadc2")){
             if ($value_found !~ /$value_test/){
               exit_error(" Method is $methodname\nbut $value_test is not in \$ricc2 group in control. Check Turbomole input.",$mdle);
             }
           }
         }
       }

       # Check CIOVERLAP
       if (($vdoth == -1) or ($vdoth == 1)){
         my $mldcio = ciocheck($mld);
	 print_STDOUT("CIOVERLAP path: $mldcio");
       }

     }

     ## ACES2 ##
     #----------------------------------------------------------------------------
     if ($progname eq "aces2") {

       # Nothing to be done.

     }

     ## MOPAC ##
     #----------------------------------------------------------------------------
     if ($progname eq "mopac") {

       # Nothing to be done	  

     }

     ## EXC_MOPAC ##
     #----------------------------------------------------------------------------
     if ($progname eq "exc_mopac") {

       # Load parameters
       ($nchrom_d,$nproc_d,$genf_d) = load_defaults_tpp($prog,"","");
       $nchrom = getkeyword($fpar, "nchrom", $nchrom_d);
       $nproc  = getkeyword($fpar, "nproc", $nproc_d);
       $genf   = getkeyword($fpar, "gen_file", $genf_d);

       # Parameter list
       $listpar = " nchrom   = $nchrom\n"
                 ." nproc    = $nproc\n"
                 ." gen_file = $genf\n";

       # Display
       if (($type eq "ic") or ($type eq "dyn")){
         print_STDOUT("$listpar\n");
       }

       # Write new parameter file
       open(INP,">$progpar") or die "$mdle Cannot write $progpar.";
       print INP "&exc_inp\n";
       print INP "$listpar";
       print INP "/\n";
       close(INP);

     }

     ## ORCA ##
     #----------------------------------------------------------------------------
     if ($progname eq "orca") {

       # Nothing to be done.

     }

     ## GAUSSIAN ##
     #----------------------------------------------------------------------------
     if ($progname eq "gau") {

       # Load parameters
       ($g_vers_d,$mocoef_d,$prt_mo_d,$td_st_d,$ld_thr_d,$kind_gau_d)=load_defaults_tpp($prog,"","");
       $mocoef   = getkeyword($progpar,"mocoef", $mocoef_d);
       $prt_mo   = getkeyword($progpar,"prt_mo", $prt_mo_d);
       $td_st    = getkeyword($progpar,"td_st", $td_st_d);
       $ld_thr   = getkeyword($progpar,"ld_thr", $ld_thr_d); #linear dependence threshold for using IOP(3/59=N)
       $kind_gau = getkeyword($progpar,"kind_gau", $kind_gau_d);
       $g_vers   = getkeyword($progpar,"g_vers", $g_vers_d);

       # Check environment variable
       if (($type eq "ic") or ($type eq "dyn")){
         $gauchk="gaussian.chk";
         $gaur  = $ENV{"g$g_vers"."root"};  # Gaussian environment
         if (!defined($gaur)){exit_error("Gaussian root variable $gaur is not defined. Check it and run again.",$mdle);}
         $gauchk="gaussian.chk";
         $gaucom="gaussian.com";
         $gaurwf="gaussian.rwf";
       }

       # Parameter list
       $listpar = " mocoef   = $mocoef\n"
                 ." prt_mo   = $prt_mo\n"
                 ." td_st    = $td_st\n"
                 ." ld_thr   = $ld_thr\n"
                 ." kind_gau = $kind_gau\n"
                 ." g_vers   = $g_vers\n";

       # Display
       if (($type eq "ic") or ($type eq "dyn")){
         print_STDOUT("$job with GAUSSIAN G$g_vers\n");
         print_STDOUT("GAUSSIAN path: $gaur\n");
         print_STDOUT("$listpar\n");
       }

       # Write new parameter file
       open(INP,">$progpar") or die "$mdle Cannot write $progpar.";
       print INP "$listpar";
       close(INP);

       # Check CIOVERLAP
       if (($vdoth == -1) or ($vdoth == 1)){
         my $mldcio = ciocheck($mld);
	 print_STDOUT("CIOVERLAP path: $mldcio");
       }

     }

     ## TINKER ##
     #----------------------------------------------------------------------------
     if (($progname eq "tinker") or ($hybrid_tpp{tinker} == 1)) {
            
       # Nothing to be done. 

     }

     ## DFTB+ ##
     #----------------------------------------------------------------------------
     if (($progname eq "dftb+") or ($hybrid_tpp{dftbp} == 1)) {

       if ($hybrid_tpp{analytical} == 1){
	  $prog = 8.5;
          %progconf   = prog_config($prog);
          $progpar    = $progconf{parfile};
          if (-s "$BASEDIR/$progpar"){system("cp -f $BASEDIR/$progpar .");}
       }

       # # Check environment
       # if (($type eq "ic") or ($type eq "dyn")){
       #   $dftbp   = $ENV{"DFTBP"};
       #   $dptools = $ENV{"DP_TOOLS"};
       #   if (!defined($dftbp)){exit_error("\$DFTBP variable is not defined. Check it and run again.",$mdle);}
       #   if (!defined($dptools)){exit_error("\$DP_TOOLS variable is not defined. Check it and run again.",$mdle);}
       # }
       $dftbp="";

       # Load parameters
       ($dftbp_exec_d,$skf_d)=load_defaults_tpp($prog,$dftbp,"");
       $dftbp_exec = getkeyword("$progpar","dftbp_exec",$dftbp_exec_d);
       $dftbp_skf  = getkeyword("$progpar","dftbp_skf",$skf_d);
       if ($dftbp_skf eq "default"){
         $dftbp_skf = $skf_d;
       }

       # Parameter list
       $listpar = " dftbp_exec = $dftbp_exec\n"
                 ." dftbp_skf  = $dftbp_skf\n";   
    
       # Display
       if (($type eq "ic") or ($type eq "dyn")){
	#print_STDOUT("DFTB+ path = $dftbp\n");
        #print_STDOUT("DP_TOOLS path = $dptools\n");
         print_STDOUT("$listpar\n");
       }

       # Write new parameter file
       open(INP,">$progpar") or die "$mdle Cannot write $progpar.";
       print INP "$listpar";
       close(INP);

       # Check CIOVERLAP
       if (($vdoth == -1) or ($vdoth == 1)){
         my $mldcio = ciocheck($mld);
	 print_STDOUT("CIOVERLAP path: $mldcio");
       }

     }

     ## DFT/MRCI
     #----------------------------------------------------------------------------
     if ($progname eq "dftci"){

       if ($type eq "dyn"){
          exit_error("DFTCI cannot be used for dynamics. Check it and run again.",$mdle); 
       }

       # Load parameter
       ($maxiter_d)=load_defaults_tpp($prog,"","");
       $maxiter = getkeyword($progpar, "maxiter", $maxiter_d);  

       # Parameter list
       $listpar = "maxiter = $maxiter\n";

       # Display
       if (($type eq "ic") or ($type eq "dyn")){
         print_STDOUT("$listpar\n");
       }

       # Write new parameter file
       open(INP,">$progpar") or die "$mdle Cannot write $progpar.";
       print INP "$listpar";
       close(INP);        
    
     }

     ## GAMESS ##
     #----------------------------------------------------------------------------
     if ($progname eq "gamess") {

       # Check environment
       if (($type eq "ic") or ($type eq "dyn")){
         $GAMESS = $ENV{"GAMESS"};
         if (!defined($GAMESS)){exit_error("\$GAMESS variable is not defined. Check it and run again.",$mdle);}
       }

       # Load parameter
       ($verno_d,$ncpus_d,$run_gamess_d,$scr_d,$mocoef_d)=load_defaults_tpp($prog,"","");
       $verno      = getkeyword("$progpar","verno",$verno_d);
       $ncpus      = getkeyword("$progpar","ncpus",$ncpus_d);
       $run_gamess = getkeyword("$progpar","run_gamess",$run_gamess_d);
       $scr        = getkeyword("$progpar","scr",$scr_d);
       $mocoef     = getkeyword("$progpar","mocoef",$mocoef_d);

       # Parameter list
       $listpar = " verno      = $verno\n"
                 ." ncpus      = $ncpus\n"
                 ." run_gamess = $run_gamess\n"
                 ." scr        = $scr\n"
                 ." mocoef     = $mocoef\n";

       # Display
       if (($type eq "ic") or ($type eq "dyn")){
         print_STDOUT("GAMESS path: $GAMESS\n");
         print_STDOUT("$listpar\n");
       }

       # Write new parameter file
       open(INP,">$progpar") or die "$mdle Cannot write $progpar.";
       print INP "$listpar";
       close(INP);

       # Define Gamess environment
       if (($type eq "ic") or ($type eq "dyn")){
         if ($scr == 0){
           $scrdir=find_scr_gamess($GAMESS,$run_gamess);
         }elsif ($scr == 1){
           $scrdir="$BASEDIR/$TP/SCR";
           make_gamess_scr("$scrdir");
         }
         print STDOUT "\nGAMESS SCR is $scrdir\n";
         test_gamess($GAMESS,$verno,$run_gamess,$scrdir);

         # modify GAMESS rungms file before each trajectory begins.
         $retval = callprogsp($mld,"create_rungms_nx.pl $run_gamess $scrdir",$mdle);
         if ($retval != 0){die "$mdle is dying now\n";}
       }

     }

     ## MLATOM ##
     #----------------------------------------------------------------------------
     if ($progname eq "mlatom") {

       # Check environment
       if (($type eq "ic") or ($type eq "dyn")){
         $mlatom  = $ENV{"MLatom"}; 
         if (!defined($mlatom)){exit_error("\$MLatom variable is not defined. Check it and run again.",$mdle);}
       }

       # Display
       if (($type eq "ic") or ($type eq "dyn")){
         print_STDOUT("\tMLATOM path: $mlatom\n");
       }

     }

     ## BAGEL ##
     #----------------------------------------------------------------------------
     if ($progname eq "bagel") {

       # Nothing to do

     }

     ## PySOC ##
     #----------------------------------------------------------------------------
     if ($progname eq "pysoc") {

       if ($lvprt < 2){exit_error("For $label, lvprt should be at least 2. Fix it and run again.",$mdle);}

     }

     ## CP2K ##
     #----------------------------------------------------------------------------
     if ($progname eq "cp2k") {

       if (($type eq "ic") or ($type eq "dyn")){

       # CP2K environment
       $cp2k  = $ENV{"CP2K"};
       $ENV{"OMP_NUM_THREADS"}=1; #for CP2K ssmp / psmp: always set OMP_NUM_THREADS=1!
                             #at least for current CP2K version

       # Load parameters
       ($parallel_d,$exec_d) = load_defaults_tpp($prog,"","");
       $parallel = getkeyword($progpar,"parallel",$parallel_d);
       $exec = getkeyword($progpar,"exec",$exec_d);

       # Check that path is exported
       if (!-e "$cp2k/$exec"){exit_error("Cannot find CP2K program at $cp2k/$exec.\nCheck \$CP2K variable.",$mdle);}

       # Parameter list
       $listpar = "parallel = $parallel\n"
                  ."exec     = $exec\n";

       # Display
       print_STDOUT("CP2K path: $cp2k\n");
       print_STDOUT("$listpar\n");

       # Write new parameter file
       open(INP,">$progpar") or die "$mdle cannot write $progpar.";
       print INP "$listpar";
       close(INP);

       # Check input
       }else{
        exit_error("CP2K can only be used for initial conditions or dynamics.",$mdle);
       }
 
       # Check CIOVERLAP
       if (($vdoth == -1) or ($vdoth == 1)){
         my $mldcio = ciocheck($mld);
	 print_STDOUT("CIOVERLAP path: $mldcio");
       }

     }

    ## FROMAGE ##
    if ($progname eq "fromage"){
      # If checkings must be done
      # include them here
    }

     ## HYBRID ##
     #----------------------------------------------------------------------------
     if ($progname eq "hybrid") {

       # Display
       if (($type eq "ic") or ($type eq "dyn")){
         print_STDOUT("$job with HYBRID gradients using:\n");
         if (-e $JAD){
           $JDIR = $JAD;
         }elsif(-e $JND){
           $JDIR = $JND;
         }else{
           exit_error("Cannot find hybrid.control file.",$mdle);
         }
         open(HC,"$JDIR/hybrid.control") or die "cannot open $JDIR/hybrid.control\n$!\n";
         while(<HC>){
           chomp; s/^\s+//;
           if(/^job /){
             @line=split/\s+/;
             print_STDOUT("\t job:\t $line[1]\n\t\t method:\t $line[3]\n\t\t regions:\t $line[2]\n\t\t factor:\t $line[4]\n");
           }
         }
       }

       # Load parameter 
       if (($type eq "ic") or ($type eq "dyn")){
         $progpar = "columbus.par";
         ($memcol_d,$cirestart_d,$reduce_tol_d,$quad_conv_d,$mc_conv_d,$ci_conv_d,$ivmode_d,$mocoef_d,$prt_mo_d,$all_grads_d,$grad_lvl_d,$cver_d,$ci_iter_d,$citol_d,$mull_pop_d)=load_defaults_tpp($prog,"","");
         $mocoef = getkeyword($progpar,"mocoef",$mocoef_d);
       }
     }

}

sub ciocheck{
#
#====================================================================================
#
# This subroutine returns $CIOVERLAP path. 
# 1- First: check $CIOVERLAP
# 2- Second: check $NX/../../CIOVERLAP 
# 3- Third: check $NX/cioverlap-64
# 4- If not found, it ends in error.
# 
#------------------------------------------------------------------------------------
#
  my $mld,$mldcio; 

  ($mld) = @_;

  if (defined $ENV{"CIOVERLAP"}){
    $mldcio = $ENV{"CIOVERLAP"};
  }elsif(-s "$mld/../../CIOVERLAP"){
    $mldcio = "$mld/../../CIOVERLAP";
  }elsif(-s "$mld/cioverlap-64"){
    $mldcio = "$mld/cioverlap-64"; 
  }else{
    exit_error("CIOVERLAP was not found. Please define \$CIOVERLAP variable.",$mdle);
  }

  return $mldcio;

}

sub read_credits{
#
#====================================================================================
#
# This subroutine returns the credits of NEWTON-X/Third-party interfaces. 
#
#------------------------------------------------------------------------------------
#
    local ($code, $value);
    ($code) = @_;
    $value = "";
    if ($code eq "gamess"){
       $value =        "-----------------------------------------------------------------------\n";
       $value = $value."NEWTON-X/GAMESS interface by A.West and T.Windus, Iowa State University\n";
       $value = $value."-----------------------------------------------------------------------\n";
    }elsif ($code eq "bagel") {
       $value=       "-----------------------------------------------------------------------\n";
       $value=$value."NEWTON-X/BAGEL interface by J.W.Park and T. Shiozaki, Northwestern University\n";
       $value=$value."-----------------------------------------------------------------------\n";
    }elsif ($code eq "fromage"){
       $value=       "--------------------------------------------------------------------------------------\n";
       $value=$value."NEWTON-X/fromage interface by Federico J. Hernandez and Rachel Crespo Otero, UCL. UK. \n";
       $value=$value."--------------------------------------------------------------------------------------\n";
    }

    return $value;
}

sub load_defaults_tpp{
#
#====================================================================================
#
# This routine contains the defaults for the third-party programs.
#
# Use it as:
# ($keyword1_d,$keyword2_d,...) = load_defaults_tpp($prog,$param1,$param2)
# where:
# $keywordi_d is the default value for keyword i.
# $prog is the value indentifying the TPP.
# $param1 and $param2 are optional parameters to control the routine behavior. 
# They depend on each TPP.
#
#------------------------------------------------------------------------------------
#
    my ($param1,$param2);
    my ($memcol_d,$cirestart_d,$reduce_tol_d,$quad_conv_d,$mc_conv_d,$ci_conv_d,$ivmode_d,$all_grads_d,$grad_lvl_d,$cver_d);
    my ($ci_iter_d,$citol_d,$mull_pop_d);
    my ($g_vers_d,$mocoef_d,$prt_mo_d,$td_st_d,$ld_thr_d,$kind_gau_d);
    my ($dftbp_exec_d,$skf_d,$dftbp);
    my ($verno_d,$ncpus_d,$run_gamess_d,$scr_d);
    my ($dftb_exec_d,$other_state_d,$mult_d);
    my ($anmod_d,$path_d,$tully_mod_d,$cs_mod_d);
    my ($k_0_d,$r1_0_d,$vr1_0_d,$k_1_d,$r1_1_d,$vr1_1_d,$k_2_d,$r1_2_d,$vr1_2_d,$vr0_1_d,$d_1_d,$alpha_1_d,$vr0_2_d,$d_2_d,$alpha_2_d,$gr_1_d,$gr_2_d,$vrc_r_d,$beta_r_d,$rc_r_d,$vrc_i_d,$beta_i_d,$rc_i_d);
    my ($A_d,$B_d,$C_d,$D_d,$Z_d,$E0_d,$theta_d,$alpha_d,$DE_d);
    my ($beta_d,$kx_d,$ky_d,$delta_d,$x1_d,$x2_d,$x3_d,$gamma_d);
    my ($e0_d,$v0_d,$N_d,$Jw_d,$sd_d);
    my ($parallel_d,$maxiter_d);
    my (%progconf,$prog,$progname,$methodname);
    my ($nchrom_d,$nproc_d,$genf_d);

    ($prog, $param1, $param2) = @_;
    %progconf   = prog_config($prog);
    $progname   = $progconf{progname};
    $methodname = $progconf{methodname};

    ## ANALYTICAL ##
    #--------------------------------------------------------------------------
    if ($progname eq "analytical"){
       if ($param1 eq "main"){   # Load defaults for the main keywords of analytical models
         $anmod_d      ="analytical.model";
         $path_d       ="\$NX";
         $tully_mod_d  = 1;
         $cs_mod_d     ="2RHE";
         return $anmod_d,$path_d,$tully_mod_d,$cs_mod_d;
       }elsif($param1 eq "cs"){  # Load defaults for CS-FSSH analytical models
         $k_0_d        = 3.0;
         $r1_0_d       = 0.0;
         $vr1_0_d      = 0.0;
         $k_1_d        = 3.0;
         $r1_1_d       = 0.25;
         $vr1_1_d      = 1.0;
         $k_2_d        = 3.0;
         $r1_2_d       = 0.25;
         $vr1_2_d      = 0.0;
         $vr0_1_d      = 2.0;
         $d_1_d        = 0.0;
         $alpha_1_d    = 0.8;
         $vr0_2_d      = 2.0;
         $d_2_d        = -1.0;
         $alpha_2_d    = 0.8;
         $gr_1_d       = 0.1;
         $gr_2_d       = 0.2;
         $vrc_r_d      = 0.03;
         $beta_r_d     = 40.0;
         $rc_r_d       = 0.464445;
         $vrc_i_d      = 0.0;
         $beta_i_d     = 0.0;
         $rc_i_d       = 0.0;
	 return $k_0_d,$r1_0_d,$vr1_0_d,$k_1_d,$r1_1_d,$vr1_1_d,$k_2_d,$r1_2_d,$vr1_2_d,$vr0_1_d,$d_1_d,$alpha_1_d,$vr0_2_d,$d_2_d,$alpha_2_d,$gr_1_d,$gr_2_d,$vrc_r_d,$beta_r_d,$rc_r_d,$vrc_i_d,$beta_i_d,$rc_i_d;
       }elsif($param1 eq "tully"){  # Load defaults for 1D collection
         if    ($param2 eq "1"){    # Model 1
	    $A_d     = 0.01;
	    $B_d     = 1.6;
	    $C_d     = 0.005;
	    $D_d     = 1.0;
	    $Z_d     = 0.0;
	    $E0_d    = 0.0;
	    $theta_d = 0.0;
	    $alpha_d = 0.0;
	    $DE_d    = 0.0;
	    return $A_d,$B_d,$C_d,$D_d,$Z_d,$E0_d,$theta_d,$alpha_d,$DE_d;
	 }elsif($param2 eq "2"){   # Model 2
	    $A_d     = 0.10;
	    $B_d     = 0.28;
	    $C_d     = 0.015;
	    $D_d     = 0.06;
	    $Z_d     = 0.05;
	    $E0_d    = 0.0;
	    $theta_d = 0.0;
	    $alpha_d = 0.0;
	    $DE_d    = 0.0;
	    return $A_d,$B_d,$C_d,$D_d,$Z_d,$E0_d,$theta_d,$alpha_d,$DE_d;
	 }elsif($param2 eq "3"){  # Model 3
	    $A_d     = 6E-4;
	    $B_d     = 0.1;
	    $C_d     = 0.9;
	    $D_d     = 0.0;
	    $Z_d     = 0.0;
	    $E0_d    = 0.0;
	    $theta_d = 0.0;
	    $alpha_d = 0.0;
	    $DE_d    = 0.0;
	    return $A_d,$B_d,$C_d,$D_d,$Z_d,$E0_d,$theta_d,$alpha_d,$DE_d;
         }elsif($param2 eq "4"){  # Model 4
	    $A_d     = 6E-4;
	    $B_d     = 0.1;
	    $C_d     = 0.9;
	    $D_d     = 0.0;
	    $Z_d     = 0.4;
	    $E0_d    = 0.0;
	    $theta_d = 0.0;
	    $alpha_d = 0.0;
	    $DE_d    = 0.0;
	    return $A_d,$B_d,$C_d,$D_d,$Z_d,$E0_d,$theta_d,$alpha_d,$DE_d;
         }elsif($param2 eq "5"){  # Model 5
	    $A_d     = 0.05;
	    $B_d     = 0.1;
	    $C_d     = 0.0;
	    $D_d     = 0.0;
	    $Z_d     = 0.0;
	    $E0_d    = 0.0;
	    $theta_d = 12;
	    $alpha_d = 2;
	    $DE_d    = 0.01;
	    return $A_d,$B_d,$C_d,$D_d,$Z_d,$E0_d,$theta_d,$alpha_d,$DE_d;
	 }
       }elsif($param1 eq "sbh"){  # Load defaults for Spin-Boson Hamiltonian
	    $e0_d = 12000;
	    $v0_d = 800;
	    $N_d  = 2;
	    $Jw_d = "user";
	    $sd_d = "138.2378     0.012\n203.879      0.000";
	    return $e0_d,$v0_d,$N_d,$Jw_d,$sd_d;   
       }elsif($param1 eq "2DCI"){ # Load defaults for 2D conical intersection
	    $alpha_d = 3.0;
            $beta_d  = 1.5;
            $kx_d    = 0.02;
            $ky_d    = 0.10;
            $delta_d = 0.01;
            $x1_d    = 4.0;
            $x2_d    = 3.0;
            $x3_d    = 3.0;
            $gamma_d = 0.04;
	    return $alpha_d,$beta_d,$kx_d,$ky_d,$delta_d,$x1_d,$x2_d,$x3_d,$gamma_d;
       }
    }
    
    ## COLUMBUS ##
    #--------------------------------------------------------------------------
    if ($progname eq "columbus"){
       $memcol_d     = 1600;
       $cirestart_d  = 0;
       $reduce_tol_d = 1;
       $quad_conv_d  = 60;
       $mc_conv_d    = 0;
       $ci_conv_d    = 0;
       $ivmode_d     = 8;
       $mocoef_d     = 1; 
       $prt_mo_d     = 0;
       $all_grads_d  = 1;
       if    ($methodname eq "columbus_mrci"){
         $grad_lvl_d   = "CI";
       }elsif($methodname eq "columbus_sa-mcscf"){
         $grad_lvl_d   = "MCSCF";
       }
       $cver_d       = 7.0;
       $ci_iter_d    = 30;
       $citol_d      = 1E-4;
       $mull_pop_d   = 0;
       return $memcol_d,$cirestart_d,$reduce_tol_d,$quad_conv_d,$mc_conv_d,$ci_conv_d,$ivmode_d,$mocoef_d,$prt_mo_d,$all_grads_d,$grad_lvl_d,$cver_d,$ci_iter_d,$citol_d,$mull_pop_d;
    }

    ## TURBOMOLE ##
    #--------------------------------------------------------------------------
    if ($progname eq "turbomole"){
       $parallel_d = 1;
       return $parallel_d;
    }

    ## EXASH MOPAC
    #--------------------------------------------------------------------------
    if ($progname eq "exc_mopac"){
      $nchrom_d = 2;
      $nproc_d = $nchrom_d;
      $genf_d = 1;
      return $nchrom_d,$nproc_d,$genf_d;
    }

    ## GAUSSIAN ##
    #--------------------------------------------------------------------------
    if ($progname eq "gau"){
       $mocoef_d   = 0; 
       $prt_mo_d   = 20;
       $td_st_d    = 0; 
       $ld_thr_d   = 14; 
       $kind_gau_d = 0;
       $g_vers_d   = "16";
       if     ($methodname eq "gaussian-cas"){
         $g_vers_d   = "09";
       }elsif ($methodname eq "gaussian-lr"){
         $g_vers_d   = "16";
       }
       return $g_vers_d,$mocoef_d,$prt_mo_d,$td_st_d,$ld_thr_d,$kind_gau_d;
    }

    ## DFTB+ ##
    #--------------------------------------------------------------------------
    if ($progname eq "dftb+"){
       $dftbp=$param1;  # param1 gets $DFTB value.
       $dftbp_exec_d="dftb+";
       $skf_d="$dftbp/sk/3ob-3-1";
       return $dftbp_exec_d,$skf_d;
    }
    
    ## DFTCI ##
    #--------------------------------------------------------------------------
    if ($progname eq "dftci"){
       $maxiter_d = 6;
       return $maxiter_d;
    }
    
    ## GAMESS ##
    #--------------------------------------------------------------------------
    if ($progname eq "gamess"){
       $verno_d      = "00";
       $ncpus_d      = 1;
       $run_gamess_d = "rungms";
       $scr_d        = 1;
       $mocoef_d     = 1;
       return $verno_d,$ncpus_d,$run_gamess_d,$scr_d,$mocoef_d;
    }

    ## CP2K ##
    #--------------------------------------------------------------------------
    if ($progname eq "cp2k"){
       $parallel_d = 1;
       $exec_d     = "cp2k.ssmp";
       return $parallel_d,$exec_d;
    }

}

sub test_list{
#
#====================================================================================
#
# List of NX tests
# 
#------------------------------------------------------------------------------------
#
    my ($all,$c_t,$turb,$col,$gau,$dftb,$anlt,$tul1,$sbh,$csa,$nx,$null,$qmmm,$td);
    my ($cc2,$adc2,$cas,$mcs,$mrc,$df,$pp,$ind,$progname,$method,$description,$dir);
    my ($games,$gamescc,$tda,$utd,$bagel,$pt2,$IS,$ST,$mlatom,$ml,$mopac,$fomo,$gaupys);
    my ($prog,$file);
    my ($reg_exp,$field,$eps,@string,$counter,$value);
    my ($e0,$e1,$e2);
    my ($m_t,$exash,$orca);
    ($ind)=@_;

    $all   ="All tests       ";
    $anlt  ="ANALYTICAL      ";
    $tul1  ="TULLY 1         ";
    $sbh   ="SPIN-BOSON      ";
    $csa   ="ANALYTIC CS-FSSH";
    $col   ="COLUMBUS        ";
    $turb  ="TURBOMOLE       ";
    $gau   ="GAUSSIAN        ";
    $gaupys="GAUSSIAN/PYSOC  ";
    $mopac ="MOPAC           ";
    $dftb  ="DFTB            ";
    $dftbp ="DFTB+           ";
    $games ="GAMESS          ";
    $bagel ="BAGEL           ";
    $mlatom="MLATOM          ";
    $c_t   ="COLUMBUS/TINKER ";
    $m_t   ="MOPAC/TINKER    ";
    $nx    ="NEWTON-X        ";
    $cp2k  ="CP2K            ";
    $orca  ="ORCA            ";

    $null    ="           ";
    $qmmm    ="QM/MM      ";
    $td      ="TDDFT      ";
    $utd     ="TD-UDFT    ";
    $tda     ="TDA        ";
    $cc2     ="RICC2      ";
    $adc2    ="ADC(2)     ";
    $cas     ="CASSCF     ";
    $mcs     ="MCSCF      ";
    $mrc     ="MRCI       ";
    $pt2     ="CASPT2     ";
    $df      ="DFTB       ";
    $ml      ="ML MODEL   ";
    $pp      ="PICK-POINTS";
    $gamescc ="CCSD       ";
    $IS      ="IMP_SAMP   ";
    $ST      ="STATISTICS ";
    $fomo    ="FOMO-CI    ";
    $exash   ="EXASH      ";

    $string[0]="%  ";                    # look field 4, eps =1E-5
    $string[1]=" Vertical excitation";   # look field 3, eps =1E-4
    $string[2]="Wave function state";    # look field 5, eps =1E-3
    $string[3]=" IP1  =";                # look field 6, eps =1E-3
    $string[4]=" 1            1 ";       # look field 4, eps =1E-3
    $string[5]="Effective number of points"; # look field 5, eps =1E-5
    $string[6]="    2.00";               # look field 3, eps =1E-5

    $e0=1E-5;
    $e1=1E-4;
    $e2=1E-3;

    # ALL TESTS
    $counter=0;
    if ($ind == $counter){
         $progname    = $all;
         $method      = $null;
         $description = "";
         $dir         = "";
         $prog        = "";
         $file        = "";
         $reg_exp     = "";
         $field       = "";
         $eps         = "";
    }

     # TURBOMOLE TDDFT AD
    $counter++;
    if ($ind == $counter){
        $progname    = $turb;
        $method      = $td;
        $description = "Dynamics, Adiabatic";
        $dir         = "MD-TM-TDDFT-AD";
        $prog        = "$mld/moldyn.pl > moldyn.log";
        $file        = "RESULTS/dyn.out";
        $reg_exp     = $string[0];
        $field       = 4;
        $eps         = $e0;
    }

    # TURBOMOLE TDDFT AND THERMOSTAT
    $counter++;
    if ($ind == $counter){
        $progname    = $turb;
        $method      = $td;
        $description = "Dynamics, Adiabatic with Andersen thermostate";
        $dir         = "MD-TM-TDDFT-THERM-AD";
        $prog        = "$mld/moldyn.pl > moldyn.log";
        $file        = "RESULTS/dyn.out";
        $reg_exp     = $string[0];
        $field       = 4;
        $eps         = $e0;
    }

    # TURBOMOLE TDDFT NAD
    $counter++;
    if ($ind == $counter){
        $progname    = $turb;
        $method      = $td;
        $description = "Dynamics, Non-adiabatic with cioverlap";
        $dir         = "MD-TM-TDDFT-NAD-CIO";
        $prog        = "$mld/moldyn.pl > moldyn.log";
        $file        = "RESULTS/dyn.out";
        $reg_exp     = $string[2];
        $field       = 5;
        $eps         = $e2;
    }

    # TURBOMOLE TDDFT NAD OD
    $counter++;
    if ($ind == $counter){
        $progname    = $turb;
        $method      = $td;
        $description = "Dynamics, Non-adiabatic with cioverlap-od";
        $dir         = "MD-TM-TDDFT-NAD-OD";
        $prog        = "$mld/moldyn.pl > moldyn.log";
        $file        = "RESULTS/dyn.out";
        $reg_exp     = $string[2];
        $field       = 5;
        $eps         = $e2;
    }

   # TURBOMOLE TDDFT
   $counter++;
   if ($ind == $counter){
       $progname    = $turb;
       $method      = $td;
       $description = "Initial conditions";
       $dir         = "IC-TM-TDDFT";
       $prog        = "$mld/initcond.pl > moldyn.log";
       $file        = "initial_condition.1.2";
       $reg_exp     = $string[1];
       $field       = 3;
       $eps         = $e1;
   }

   # TURBOMOLE RICC2
   $counter++;
   if ($ind == $counter){
       $progname    = $turb;
       $method      = $cc2;
       $description = "Dynamics, Adiabatic";
       $dir         = "MD-TM-RICC2-AD";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[0];
       $field       = 4;
       $eps         = $e0;
   }

   # TURBOMOLE RICC2
   $counter++;
   if ($ind == $counter){
       $progname    = $turb;
       $method      = $cc2;
       $description = "Dynamics, Non-adiabatic with cioverlap";
       $dir         = "MD-TM-RICC2-NAD-CIO";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # TURBOMOLE RICC2
   $counter++;
   if ($ind == $counter){
       $progname    = $turb;
       $method      = $cc2;
       $description = "Initial conditions";
       $dir         = "IC-TM-RICC2";
       $prog        = "$mld/initcond.pl > moldyn.log";
       $file        = "initial_condition.1.2";
       $reg_exp     = $string[1];
       $field       = 3;
       $eps         = $e1;
   }

   # TURBOMOLE ADC(2) NAD CIO
   $counter++;
   if ($ind == $counter){
       $progname    = $turb;
       $method      = $adc2;
       $description = "Dynamics, Non-adiabatic with cioverlap";
       $dir         = "MD-TM-ADC2-NAD-CIO";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # TURBOMOLE MP2 LP-ZPE
   $counter++;
   if ($ind == $counter){
       $progname    = $turb;
       $method      = "MP2        ";
       $description = "Dynamics, adiabatic with LP-ZPE";
       $dir         = "MD-TM-MP2-AD-LPZPE";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # INITIAL CONDITIONS WITH SAMP_DIST
   $counter++;
   if ($ind == $counter){
       $progname    = $turb;
       $method      = $adc2;
       $description = "Initial conditions with distribution transformations";
       $dir         = "IC-TM-ADC2-DISTRIB";
       $prog        = "$mld/initcond.pl > moldyn.log";
       $file        = "initial_condition.1.3";
       $reg_exp     = $string[1];
       $field       = 3;
       $eps         = $e1;
   }

   # TURBOMOLE ADC(2) NAD OD
   $counter++;
   if ($ind == $counter){
       $progname    = $turb;
       $method      = $adc2;
       $description = "Dynamics, Non-adiabatic with cioverlap-od";
       $dir         = "MD-TM-ADC2-NAD-OD";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # COLUMBUS CASSCF AD
   $counter++;
   if ($ind == $counter){
       $progname    = $col;
       $method      = $cas;
       $description = "Dynamics, Adiabatic";
       $dir         = "MD-COL-CASSCF-AD";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[0];
       $field       = 4;
       $eps         = $e0;
   }

   # COLUMBUS MCSCF NAD VECTORS
   $counter++;
   if ($ind == $counter){
       $progname    = $col;
       $method      = $mcs;
       $description = "Dynamics, Non-adiabatic with NAC vectors and 3 states";
       $dir         = "MD-COL-MCSCF-NAD-VEC";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # COLUMBUS MRCI NAD VECTORS
   $counter++;
   if ($ind == $counter){
       $progname    = $col;
       $method      = $mrc;
       $description = "Dynamics, Non-adiabatic with NAC vectors and 2 states";
       $dir         = "MD-COL-MRCI-NAD-VEC";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # COLUMBUS MRCI NAD CIOVERLAP
   $counter++;
   if ($ind == $counter){
       $progname    = $col;
       $method      = $mrc;
       $description = "Dynamics, Non-adiabatic with ciovelap";
       $dir         = "MD-COL-MRCI-NAD-CIO";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # COLUMBUS MRCI NAD LOCAL-DIABATIZATION
   $counter++;
   if ($ind == $counter){
       $progname    = $col;
       $method      = $mrc;
       $description = "Dynamics, Non-adiabatic with local-diabatization";
       $dir         = "MD-COL-MRCI-NAD-LD";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # COLUMBUS MCSCF NAD TD-BA
   $counter++;
   if ($ind == $counter){
       $progname    = $col;
       $method      = $cas;
       $description = "Dynamics, Non-adiabatic with TD-BA";
       $dir         = "MD-COL-MCSCF-NAD-TDBA";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # COLUMBUS MRCI CS-FSSH MODEL 1
   $counter++;
   if ($ind == $counter){
       $progname    = $col;
       $method      = $mrc;
       $description = "Dynamics, CS-FSSH gamma model 1";
       $dir         = "MD-COL-MRCI-NAD-VEC-CSFSSH-1";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # COLUMBUS MRCI CS-FSSH MODEL 2
   $counter++;
   if ($ind == $counter){
       $progname    = $col;
       $method      = $mrc;
       $description = "Dynamics, CS-FSSH gamma model 2";
       $dir         = "MD-COL-MRCI-NAD-VEC-CSFSSH-2";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # COLUMBUS/TINKER QM/MM
   $counter++;
   if ($ind == $counter){
       $progname    = $c_t;
       $method      = $qmmm;
       $description = "Dynamics, Non-adiabatic with QM/MM";
       $dir         = "MD-COL-TNK-QMMM-NAD-VEC";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # COLUMBUS MCSCF
   $counter++;
   if ($ind == $counter){
       $progname    = $col;
       $method      = $mcs;
       $description = "Initial conditions MCSCF";
       $dir         = "IC-COL-MCSCF";
       $prog        = "$mld/initcond.pl > moldyn.log";
       $file        = "initial_condition.1.2";
       $reg_exp     = $string[1];
       $field       = 3;
       $eps         = $e1;
   }

   # COLUMBUS MRCI
   $counter++;
   if ($ind == $counter){
       $progname    = $col;
       $method      = $cas;
       $description = "Initial conditions MRCI";
       $dir         = "IC-COL-MRCI";
       $prog        = "$mld/initcond.pl > moldyn.log";
       $file        = "initial_condition.1.2";
       $reg_exp     = $string[1];
       $field       = 3;
       $eps         = $e1;
   }

   # DFTB-0
   $counter++;
   if ($ind == $counter){
       $progname    = $dftb;
       $method      = $df;
       $description = "Dynamics, Adiabatic in excited state";
       $dir         = "MD-DFTB0-AD";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[0];
       $field       = 4;
       $eps         = $e0;
   }

   # DFTB+
   $counter++;
   if ($ind == $counter){
       $progname    = $dftbp;
       $method      = $df;
       $description = "Dynamics, Nonadiabatic with cioverlap_od";
       $dir         = "MD-DFTB+-TDDFTB-NAD-OD";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # DFTB+
   $counter++;
   if ($ind == $counter){
       $progname    = $dftbp;
       $method      = $df;
       $description = "Initial conditions";
       $dir         = "IC-DFTB+-TDDFTB";
       $prog        = "$mld/initcond.pl > moldyn.log";
       $file        = "initial_condition.1.2";
       $reg_exp     = $string[1];
       $field       = 3;
       $eps         = $e1;
   }

   # GAUSSIAN CASSCF NAD
   $counter++;
   if ($ind == $counter){
       $progname    = $gau;
       $method      = $cas;
       $description = "Dynamics, Non-adiabatic with NAC vectors";
       $dir         = "MD-GAU-CASSCF-NAD-VEC";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # GAUSSIAN TDDFT NAD
   $counter++;
   if ($ind == $counter){
       $progname    = $gau;
       $method      = $td;
       $description = "Dynamics, Non-adiabatic with cioverlap";
       $dir         = "MD-GAU-TDDFT-NAD-CIO";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # GAUSSIAN TDA NAD
   $counter++;
   if ($ind == $counter){
       $progname    = $gau;
       $method      = $tda;
       $description = "Dynamics, Non-adiabatic with cioverlap";
       $dir         = "MD-GAU-TDA-NAD-CIO";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # GAUSSIAN TD-UDFT NAD
   $counter++;
   if ($ind == $counter){
       $progname    = $gau;
       $method      = $utd;
       $description = "Dynamics, Non-adiabatic with cioverlap";
       $dir         = "MD-GAU-TDUDFT-NAD-CIO";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # GAUSSIAN TDDFT NAD-OD
   $counter++;
   if ($ind == $counter){
       $progname    = $gau;
       $method      = $td;
       $description = "Dynamics, Non-adiabatic with cioverlap_od";
       $dir         = "MD-GAU-TDDFT-NAD-OD";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # GAUSSIAN TDDFT INITIAL CONDITIONS
   $counter++;
   if ($ind == $counter){
       $progname    = $gau;
       $method      = $td;
       $description = "Initial conditions";
       $dir         = "IC-GAU-TDDFT";
       $prog        = "$mld/initcond.pl > moldyn.log";
       $file        = "initial_condition.1.2";
       $reg_exp     = $string[1];
       $field       = 3;
       $eps         = $e1;
   }

   # GAUSSIAN TDDFT/PYSOC INITIAL CONDITIONS
   $counter++;
   if ($ind == $counter){
       $progname    = $gaupys;
       $method      = $td;
       $description = "Initial conditions with spin-orbit couplings";
       $dir         = "IC-GAU-PYSOC-TDDFT";
       $prog        = "$mld/initcond.pl > moldyn.log";
       $file        = "initial_condition.1.2";
       $reg_exp     = $string[1];
       $field       = 3;
       $eps         = $e1;
   }

   # GAUSSIAN TDDFT DYSON ORBITALS
   $counter++;
   if ($ind == $counter){
       $progname    = $gau;
       $method      = $td;
       $description = "Dyson orbitals";
       $dir         = "IC-GAU-TDDFT-DO";
       $prog        = "$mld/initcond.pl > moldyn.log";
       $file        = "ip.dat";
       $reg_exp     = $string[3];
       $field       = 6;
       $eps         = $e2;
   }

   # GAUSSIAN TDDFT PE SPECTRUM
   $counter++;
   if ($ind == $counter){
       $progname    = $gau;
       $method      = $td;
       $description = "Photoelectron Spectrum";
       $dir         = "SPECTRUM-GAU-TDDFT-PE";
       $prog        = "$mld/ixsec.pl > moldyn.log";
       $file        = "general.dat";
       $reg_exp     = $string[4];
       $field       = 4;
       $eps         = $e2;
   }

   # GAMESS CASSCF INITIAL CONDITIONS
   $counter++;
   if ($ind == $counter){
       $progname    = $games;
       $method      = $cas;
       $description = "Initial conditions";
       $dir         = "IC-GAMESS-CASSCF";
       $prog        = "$mld/initcond.pl > moldyn.log";
       $file        = "initial_condition.1.2";
       $reg_exp     = $string[1];
       $field       = 3;
       $eps         = $e1;
   }

   # GAMESS CASSCF NAD
   $counter++;
   if ($ind == $counter){
       $progname    = $games;
       $method      = $cas;
       $description = "Dynamics, Non-adiabatic with NAC vectors";
       $dir         = "MD-GAMESS-CASSCF-NAD-VEC";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # GAMESS CASSCF AD
   $counter++;
   if ($ind == $counter){
       $progname    = $games;
       $method      = $cas;
       $description = "Dynamics, Adiabatic";
       $dir         = "MD-GAMESS-CASSCF-AD";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[0];
       $field       = 4;
       $eps         = $e0;
   }

   # GAMESS CCSD AD
   $counter++;
   if ($ind == $counter){
       $progname    = $games;
       $method      = $gamescc;
       $description = "Dynamics, Adiabatic in excited state";
       $dir         = "MD-GAMESS-CCSD-AD";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[0];
       $field       = 4;
       $eps         = $e0;
   }

   # MOPAC FOMO-CI INITIAL CONDITIONS
   $counter++;
   if ($ind == $counter){
       $progname    = $mopac;
       $method      = $fomo;
       $description = "Initial conditions";
       $dir         = "IC-MOPAC-FOMOCI";
       $prog        = "$mld/initcond.pl > moldyn.log";
       $file        = "initial_condition.1.2";
       $reg_exp     = $string[1];
       $field       = 3;
       $eps         = $e1;
   }   

   # MOPAC FOMO-CI NAD LD
   $counter++;
   if ($ind == $counter){
       $progname    = $mopac;
       $method      = $fomo;
       $description = "Dynamics, nonadiabatic with local diabatization";
       $dir         = "MD-MOPAC-FOMOCI-NAD-LD";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[0];
       $field       = 4;
       $eps         = $e0;
   }

   # MOPAC FOMO-CI NAD CIO
   $counter++;
   if ($ind == $counter){
       $progname    = $mopac;
       $method      = $fomo;
       $description = "Dynamics, nonadiabatic with time derivative couplings";
       $dir         = "MD-MOPAC-FOMOCI-NAD-CIO";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[0];
       $field       = 4;
       $eps         = $e0;
   }

   # MOPAC FOMO-CI NAD VEC
   $counter++;
   if ($ind == $counter){
       $progname    = $mopac;
       $method      = $fomo;
       $description = "Dynamics, nonadiabatic with nonadiabatic vector";
       $dir         = "MD-MOPAC-FOMOCI-NAD-VEC";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[0];
       $field       = 4;
       $eps         = $e0;
   }

   # MOPAC FOMO-CI NAD TDBA
   $counter++;
   if ($ind == $counter){
       $progname    = $mopac;
       $method      = $fomo;
       $description = "Dynamics, nonadiabatic with TD-BA";
       $dir         = "MD-MOPAC-FOMOCI-NAD-TDBA";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[0];
       $field       = 4;
       $eps         = $e0;
   }

   # MOPAC/TINKER EXASH NAD
   $counter++;
   if ($ind == $counter){
       $progname    = $m_t;
       $method      = $exash;
       $description = "Dynamics, nonadiabatic with local diabatization";
       $dir         = "MD-MOPAC-TNK-EXASH-NAD-LD";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[0];
       $field       = 4;
       $eps         = $e0;
   }

   # BAGEL CASSCF INITIAL CONDITIONS
   $counter++;
   if ($ind == $counter){
       $progname    = $bagel;
       $method      = $cas;
       $description = "Initial conditions";
       $dir         = "IC-BAGEL-CASSCF";
       $prog        = "$mld/initcond.pl > moldyn.log";
       $file        = "initial_condition.1.2";
       $reg_exp     = $string[1];
       $field       = 3;
       $eps         = $e1;
   }

   # BAGEL CASSCF NAD
   $counter++;
   if ($ind == $counter){
       $progname    = $bagel;
       $method      = $cas;
       $description = "Dynamics, Non-adiabatic with NAC vectors";
       $dir         = "MD-BAGEL-CASSCF-NAD-VEC";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # BAGEL CASSCF NAD
   $counter++;
   if ($ind == $counter){
       $progname    = $bagel;
       $method      = $pt2;
       $description = "Dynamics, Non-adiabatic with NAC vectors";
       $dir         = "MD-BAGEL-CASPT2-NAD-VEC";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # ML ATOM NAD WITH TDBA
   $counter++;
   if ($ind == $counter){
       $progname    = $mlatom;
       $method      = $ml;
       $description = "Dynamics, Non-adiabatic with TD-BA";
       $dir         = "MD-MLATOM-NAD-TDBA";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # CP2K TDDFT AD
    $counter++;
    if ($ind == $counter){
        $progname    = $cp2k;
        $method      = $td;
        $description = "Dynamics, Adiabatic";
        $dir         = "MD-CP2K-TDDFT-GAPW-AD";
        $prog        = "$mld/moldyn.pl > moldyn.log";
        $file        = "RESULTS/dyn.out";
        $reg_exp     = $string[0];
        $field       = 4;
        $eps         = $e0;
    }

   # CP2K TDDFT AD
    $counter++;
    if ($ind == $counter){
        $progname    = $cp2k;
        $method      = $td;
        $description = "Dynamics, Adiabatic";
        $dir         = "MD-CP2K-TDDFT-GPW-AD";
        $prog        = "$mld/moldyn.pl > moldyn.log";
        $file        = "RESULTS/dyn.out";
        $reg_exp     = $string[0];
        $field       = 4;
        $eps         = $e0;
    }

    # CP2K TDDFT AD AND THERMOSTAT
    $counter++;
    if ($ind == $counter){
        $progname    = $cp2k;
        $method      = $td;
        $description = "Dynamics, Adiabatic with Andersen thermostate";
        $dir         = "MD-CP2K-TDDFT-GAPW-THERM-AD";
        $prog        = "$mld/moldyn.pl > moldyn.log";
        $file        = "RESULTS/dyn.out";
        $reg_exp     = $string[0];
        $field       = 4;
        $eps         = $e0;
    }

    # CP2K TDDFT NAD OD
    $counter++;
    if ($ind == $counter){
        $progname    = $cp2k;
        $method      = $td;
        $description = "Dynamics, Non-adiabatic with cioverlap-od";
        $dir         = "MD-CP2K-TDDFT-NAD-OD";
        $prog        = "$mld/moldyn.pl > moldyn.log";
        $file        = "RESULTS/dyn.out";
        $reg_exp     = $string[2];
        $field       = 5;
        $eps         = $e2;
    }    

   # CP2K TDDFT Initial conditions
   $counter++;
   if ($ind == $counter){
       $progname    = $cp2k;
       $method      = $td;
       $description = "Initial conditions";
       $dir         = "IC-CP2K-TDDFT-GAPW";
       $prog        = "$mld/initcond.pl > moldyn.log";
       $file        = "initial_condition.1.2";
       $reg_exp     = $string[1];
       $field       = 3;
       $eps         = $e1;
   }

   # CP2K TDDFT Initial conditions
   $counter++;
   if ($ind == $counter){
       $progname    = $cp2k;
       $method      = $td;
       $description = "Initial conditions";
       $dir         = "IC-CP2K-TDDFT-GPW";
       $prog        = "$mld/initcond.pl > moldyn.log";
       $file        = "initial_condition.1.2";
       $reg_exp     = $string[1];
       $field       = 3;
       $eps         = $e1;
   }

   # ORCA TDDFT NAD WITH TDBA
   $counter++;
   if ($ind == $counter){
       $progname    = $orca;
       $method      = $td;
       $description = "Dynamics, Non-adiabatic with TD-BA";
       $dir         = "MD-ORCA-TDDFT-NAD-TDBA";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }   

   # ORCA TDDFT INITIAL CONDITIONS
   $counter++;
   if ($ind == $counter){
       $progname    = $orca;
       $method      = $td;
       $description = "Initial conditions";
       $dir         = "IC-ORCA-TDDFT";
       $prog        = "$mld/initcond.pl > moldyn.log";
       $file        = "initial_condition.1.2";
       $reg_exp     = $string[1];
       $field       = 3;
       $eps         = $e1;
   }

   # ANALYTICAL MODEL 1
   $counter++;
   if ($ind == $counter){
       $progname    = $nx;
       $method      = $anlt;
       $description = "Dynamics, Non-adiabatic with NAC vectors using model 1";
       $dir         = "MD-ANALYT-MODEL1-NAD-VEC";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # ANALYTICAL TULLY MODEL 1
   $counter++;
   if ($ind == $counter){
       $progname    = $nx;
       $method      = $tul1;
       $description = "Dynamics, Non-adiabatic with NAC vectors using Tully model 1";
       $dir         = "MD-ANALYT-TULLY1-NAD-VEC";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # ANALYTICAL SPIN-BOSON HAMILTONIAN
   $counter++;
   if ($ind == $counter){
       $progname    = $nx;
       $method      = $sbh;
       $description = "Dynamics, Non-adiabatic with NAC vectors using Spin-Boson Hamiltonian";
       $dir         = "MD-ANALYT-SPINBOSON-NAD-VEC";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # ANALYTICAL CS-FSSH HAMILTONIAN 2RHE
   $counter++;
   if ($ind == $counter){
       $progname    = $nx;
       $method      = $csa;
       $description = "Dynamics, CS-FSSH using 2RHE Hamiltonian";
       $dir         = "MD-ANALYT-CSFSSH-2RHE";
       $prog        = "$mld/moldyn.pl > moldyn.log";
       $file        = "RESULTS/dyn.out";
       $reg_exp     = $string[2];
       $field       = 5;
       $eps         = $e2;
   }

   # PICK POINTS
   $counter++;
   if ($ind == $counter){
       $progname    = $nx;
       $method      = $pp;
       $description = "Initial conditions from a previous dynamics";
       $dir         = "IC-PICK-POINTS";
       $prog        = "$mld/initcond.pl > moldyn.log";
       $file        = "initial_condition.3.1";
       $reg_exp     = $string[1];
       $field       = 3;
       $eps         = $e1;
   }

   # fk
   # SPECTRUM WITH IMPORTANCE SAMPLING
   $counter++;
   if ($ind == $counter){
       $progname    = $nx;
       $method      = $IS;
       $description = "Absorption spectrum with importance sampling";
       $dir         = "SPECTRUM-IMPORTANCE_SAMPLING";
       $prog        = "$mld/makedir.pl > moldyn.log";
       $file        = "importance_sampling.log";
       $reg_exp     = $string[5];
       $field       = 5;
       $eps         = $e0;
   }

   # STATISTICAL ANALYSIS
   $counter++;
   if ($ind == $counter){
       $progname    = $nx;
       $method      = $ST;
       $description = "Statistical analysis of dynamics";
       $dir         = "DYN";
       $prog        = "$mld/analysis.pl > moldyn.log";
       $file        = "ANALYSIS/mean_value.1";
       $reg_exp     = $string[6];
       $field       = 3;
       $eps         = $e0;
   }

    if ($ind == -1){
         $value=$counter;
         return $value;
    }else{
         return $progname,$method,$description,$dir,$prog,$file,$reg_exp,$field,$eps;
   }

}

sub osc_strength{
#
#====================================================================================
#
# This subroutine reads and returns oscillator strengths from TPP to NX.
# Usage $OOS=osc_strength(mdle,nis,nfs,prog)
# where mdle is the program calling the routine, nis and nfs are the initial
# and final states, and prog is the quantum chemistry program used.
#
# For new TPP interfaces, it is recommended to add rhe oos reading directly in
# the run_tpp_initcond.pl and run-tpp.pl programs, not here.
#
#------------------------------------------------------------------------------------
#
    local ($mdle,$nis,$nfs,$mxns,$prog,$value,$found,$file,$g,$dx,$dy,$dz,$f,$i,$j,$ct,$aux,@h,$si,$sf);
    local (%progconf);
    local ($hybcolumbus);

    ($mdle,$nis,$nfs,$prog)=@_;

    $nis=$nis*1;
    $nfs=$nfs*1;
    $found = 0;

    if ($nis < $nfs){
      $mxns=$nfs;
    }elsif($nis > $nfs){
      $mxns=$nis;
    }

    %progconf   = prog_config($prog); 
    $progname   = $progconf{progname};
    $methodname = $progconf{methodname};

    $hybcolumbus="n";
    if ($methodname eq "hybrid"){
       if (-e "control.run"){
          $hybcolumbus="y";
       }
    }

    if ($progname eq "columbus" or $hybcolumbus eq "y"){
    # COLUMBUS
      $file = "LISTINGS/trncils\.FROMdrt1\.state$nfs"."TOdrt1\.state$nis";
      if (!-e "$file"){
         $file = "LISTINGS/trncils\.FROMdrt1\.state$nis"."TOdrt1\.state$nfs";
      }
      # print_STDOUT("$mdle Looking for oscillator strength at $file \n");
      if (-e "$file"){
        open(FL,"$file") or die "Cannot open $file to read!";
        while(<FL>){
          if (/Transition moment components/){
            while(<FL>){
               if (/electron/){
                  chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
                  ($g,$dx,$dy,$dz)=split(/\s+/,$_);
                  while(<FL>){
                     if (/Oscillator strength :/){
                        chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
                        ($g,$g,$g,$f)=split(/\s+/,$_);
                        $value="$f,$dx,$dy,$dz";
                        $found = 1;
                        last;
                     }
                  }
                  last;  
               }
            }
            last;
          }
        }
        close(FL);
      }elsif (!-e "$file"){
         $aux=0;
         $file = "LISTINGS/transls\.FROMdrt1\.state$nfs"."TOdrt1\.state$nis";
         if (!-e "$file"){
            $file = "LISTINGS/transls\.FROMdrt1\.state$nis"."TOdrt1\.state$nfs";
            if (-e $file){
               $aux=1;
            }
         }elsif(-e "$file"){
            $aux=1;
         }
         if ($aux==1){
            open(FL,"$file") or die "Cannot open $file to read!";
            while(<FL>){
              if (/Transition moment components/){
                 while(<FL>){
                   if (/bohr/){
                      chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
                      ($dx,$dy,$dz)=split(/\s+/,$_);
                      while(<FL>){
                         if (/Oscillator strength:/){
                            chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
                            ($g,$g,$f)=split(/\s+/,$_);
                            $value="$f,$dx,$dy,$dz";
                            $found = 1;
                            last;
                         }
                      }
                      last;
                   }
                }
                last;
              }
            }
         }
      }
    }

    if (($methodname eq "turbomole-ricc2") or ($methodname eq "turbomole-riadc2")){
    # TURBOMOLE RI-CC2
       $file = "grad.out";
       open(FL,"$file") or die "Cannot open $file to read";
       while(<FL>){
          if (/oscillator strength \(length gauge\)/){
             chomp;
             $found++;
             if ($found == $mxns-1){
               ($g,$g,$g,$g,$g,$g,$value)=split(/\s+/,$_);
               last;
             }
          }
       }
       close(FL);
    }

    if ($methodname eq "turbomole-tddft"){
    # TURBOMOLE TD-DFT
       $file = "grad.out";
       open(FL,"$file") or die "Cannot open $file to read!";
       while(<FL>){
          if (/Oscillator strength/){
              while(<FL>){
                  if (/length representation/){
                     chomp;
                     $found++;
                     if ($found == 2*$mxns-3){
                       ($grb,$grb,$grb,$value)=split(/\s+/,$_);
                       last;
                     }
                  }
              }
          }
       }
       close(FL);
     }

    if ($progname eq "aces2"){
    # ACES2
       $file = "aces2ls";
       open(FL,"$file") or die "Cannot open $file to read!";
       while(<FL>){
          if (/Strength/){
             $found++;
             if ($found == $mxns-1){
               $_=<FL>;
               $_=<FL>;
               chomp;
               ($grb,$grb,$grb,$grb,$dx,$dy,$dz,$f)=split(/\s+/,$_);
               $value="$f,$dx,$dy,$dz";
               $value=~s/d/e/gi;
               last;
             }
          }
       }
       close(FL);
    }

    if ($progname eq "dftci"){
    # DFT-MRCI
       $file = "output.prp";
       open(FL,"$file") or die "Cannot open $file to read!";
       while(<FL>){
          if (/T - D E N S I T Y - O U T P U T/){
             while(<FL>){
               if (/excitation E/){
                 chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
                 @h=split(/\s+/,$_);
                 $h[1] =~ s/a//;
                 $si = $h[1];
                 $h[3] =~ s/a//;
                 $sf = $h[3]; # excitation si -> sf
                 if ((($si == $nis) and ($sf == $nfs)) or
                     (($si == $nfs) and ($sf == $nis))){
                   while(<FL>){
                      if (/dipole \(L\)/){
                        chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
                        ($g,$g,$g,$dx,$dy,$dz)=split(/\s+/,$_);  #dipole
                        while(<FL>){
                          if (/osc\.str\. \(L\)/){
                            chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
                            ($g,$g,$f)=split(/\s+/,$_);
                            $value="$f,$dx,$dy,$dz";
                            $found=1;
                            last;
                          }
                        }
                        last;
                      }
                   }
                   last;
                 }
               }
             }
             last;
          }
       }
    }

    if ($progname eq "gamess"){
       # GAMESS
       #
       $ct = 0;
       for ($i = 2; $i <= $nis; $i++){
         for ($j = 1 ; $j <= $i-1; $j++) {
           $ct++;
           last if ($j==$nfs && $i==$nis);
         }
       }
       #
       $file = "osc_dipole_store";
       open(FL,"$file") or die "Cannot open $file to read";
       while(<FL>){
         # found set to 0 at very beginning of this sub
         $found++;
         if ($found == $ct){
           chomp;
           $value = $_;
           last;
         }
       }
       close(FL);
    }

    if ($progname eq "bagel") {
       # BAGEL
       #
       $ct = 0;
       for ($i = 2; $i <= $nis; $i++){
         for ($j = 1 ; $j <= $i-1; $j++) {
           $ct++;
           last if ($j==$nfs && $i==$nis);
         }
       }
       #
       $file = "bagel.out";
       system ("grep 'Oscillator strength' $file > osc_dipole_store");
       $file = "osc_dipole_store";
       open(FL,"$file") or die "We do not have $file to read oscillator strength";
       while(<FL>) {
         $found++;
         if ($found == $ct) {
           chomp;
           $_ =~ s/^\s*//;         # remove leading blanks
           $_ =~ s/\s*$//;         # remove trailing blanks
           ($g,$g,$g,$g,$g,$g,$g,$g,$g,$value)=split(/\s+/,$_);
           last;
         }
       }
       close(FL);
    }

    if ($found == 0){
    # If f is not found, return f = 1.
      print_STDOUT("\n$mdle Oscillator strength ($nis -> $nfs) was not found.
It will be assumed to be 1.0.\n");
      $value = 1.0;
    }

    return $value;
}

#
#====================================================================================
#       COLUMBUS-RELATED ROUTINES
#====================================================================================
#

sub columbus_mode{
#
#====================================================================================
#
# Look at control.run and return
#   mode = ss if ciudg
#   mode = ms if ciudgav
#   mode = mc if no ci program is found
#
#------------------------------------------------------------------------------------
#
  my ($value,$found1,$found2);
  $found1=searchkeyword("control.run","ciudgav")+searchkeyword("control.run","pciudgav");
  if ($found1 == 0){
     $found2=searchkeyword("control.run","ciudg")+searchkeyword("control.run","pciudg");
     if ($found2 == 0){
        $found3=searchkeyword("control.run","mcscfgrad");
        if ($found3 == 0){
           $value = "mc";
        }
        else
        {
           $value = "mcs";
        }
     }elsif($found2 > 0){
        $value = "ss";
     }
  }elsif($found1 > 0){
     $value = "ms";
  }
  return $value;
}

sub columbus_mem{
#
#====================================================================================
#
# Return memory value for differentColumbus versions.
#
#------------------------------------------------------------------------------------
#   
  my ($value,$cver,$memcol);

  ($cver,$memcol)=@_;

  if     ( $cver < 6.0 ){                # in 5.9.2, Columbus core memory is given in words
     $value = $memcol/7.63E-6;
  }elsif (($cver >= 6) and ($cver < 7)){ # in version 6, mem is given in Mwords
     $value = $memcol/7.63;
  }elsif ($cver >= 7) {                  # in version 7+, mem is given in MB
     $value = $memcol;
  }
  return $value;

}

sub lookatcolumbus{
#
#====================================================================================
#
# Call test of Columbus version.
#
#------------------------------------------------------------------------------------
#
          $file = $_[0];
          open(RE,"$file") or warn "Cannot open $file!"; # I'm not sure: die or warn?
          while(<RE>){
             if (/bummer/){
               if (!/warning/i){
                 print_STDOUT("ERROR: bummer found in $file. \nThe message is: $_\n");
                 $status = "with ERROR";
                 test_columbus_version();
               }
             }
          }
          close(RE);
}

sub mcscfconv{
#
#====================================================================================
#
# Return convergence infor for MCSCF.
#
#------------------------------------------------------------------------------------
#

  ($filein)=@_;
  my ($dum1,$dum2);
    $nconv=0;
      open(FILE,"$filein") or die" cannot open file: $filein\n";
      $/="\n";
         while (<FILE>)
         {
            if (/final mcscf/)
             {
              $_=<FILE> ; chop ; s/^ *//;
              ($dum1,$dum2)=split(/\s+/,$_,3);
              if ($dum1 eq "iter="){$niter=$dum2;}
              else {$niter=$dum1;}
              $nconv = /not conv/;
             }
         }
    close FILE;
    return $nconv,$niter;
} 

sub ciconv{
#
#====================================================================================
#
# Return convergence info for MRCI.
#
#------------------------------------------------------------------------------------
#

  ($filename)=@_;
  my ($dum,$niter,$nconv);
    $nconv=0;
      open(FILE,"$filename") or die" cannot open file: $filename\n";
      $/="\n";
         while (<FILE>)
         {
            if (/iterations/)
             {
              chop ; s/^ *//;
              (@dum[0..4],$niter)=split(/\s+/,$_,7);
              $nconv = /not reached/;
             }
         }
    close FILE;
    return $nconv,$niter;
}

sub test_columbus_version{
#
#====================================================================================
#
# Return Columbus version. 
#
#------------------------------------------------------------------------------------
#
  my ($nstat,$nstatdyn,$dum,$cbus,$g,$mode,$ctd,@status);
  $cbus   = $ENV{"COLUMBUS"};
  $mode   = columbus_mode();
  $ctd    = "control.d";
  if ($mode eq "ss"){
    @status   = load_status($ctd,$mdle);
    $nstat    = $status[2];
    $nstatdyn = $status[3];
    if ($nstat > $nstatdyn){
       $g=qx(grep -c "Gradient calculated for root" $cbus/runc);
       chomp($g);$g =~ s/^\s*//;$g =~ s/\s*$//;
       if ($g ne "1"){
         open(IE,">inderr");
         print IE "671";
         close(IE);
       }
    }
  }
}

sub frootci {
#
#====================================================================================
#
# Extracts the root followed in CI from ciudgls*
#
#------------------------------------------------------------------------------------
#
     local ($energy,$found);
     ($filename)=@_;

      open (ANYFILE, $filename) or die " Error in open file $filename\n";
      $/="\n";
       while ( <ANYFILE> )
        {if ( /vector at position/ )
          { $found=$_; last;}
        }
        chop $found;
       $energy=$found;
       $energy=~s/^.*energy=//g;
       $energy=~s/ //g;
       $found =~ s/energy.*$//g;
       $found =~ s/[ a-zA-Z]//g;
       return ($found,$energy);
}

sub getkeydrt {
#
#====================================================================================
#
# Get the value of a variable in cidrt and mcdrt files of COLUMBUS.
# Usage: $keyword=getkeyword(<file>,<pattern>,<default>);
# In these files, the variable appears in the first field, e.g.,
#  0 / input the maximum excitation level
# This subroutine applied to the previous line in ciudgin file:
# $mel = getkeydrt("cidrtin","maximum excitation level","2");
# returns
# $mel = 0
# If the pattern is not found or the file does not exist, value
# assumes the default value.
#
#------------------------------------------------------------------------------------
#
    local ($file,$pattern,$default,$found,$value,$g);
    ($file,$pattern,$default)=@_;
    $found = 0;
    if (-s $file){
       open(FR,"$file") or die "Cannot open $file to read!";
       while(<FR>){
          if (/\b$pattern\b/i){
            ($value,$g)=split(/\//,$_);
            $value =~ s/^\s*//;         # remove leading blanks
            $value =~ s/\s*$//;         # remove trailing blanks
            $found = 1;
          }
       }
       close(FR);
    }
    if ($found == 0){
       $value = $default;
    }
    return $value;
  }

sub wf_info{
#
#====================================================================================
#
# Get wavefunction information in third party programs.
# Only implemented for Columbus.
# Usage:
# ($c,$c2,$vec,$more,$nelem)=wf_info($state,$max_conf,$file,$progname);
#   state    = state of interest
#   max_conf = max. number of configurations to be searched for
#   file     = output file containing the information
#   progname = third party program as given by &progname($prog)
# Example: Most important (first) cfs in state 2.
# ($c,$c2,$vec,$more,$nelem)=wf_info(2,1,"LISTINGS/cipcls","columbus");
#   c        = coefficient
#   c2       = squared coeff.
#   vec      = step vector
#   more     = other information (type of csf for Columbus)
#   nelem    = number of elements found (may be smaller than max_conf)
# The routine returns the reference to arrays. The elements may be accessed
# for example:
# print "csf 1 c = $$c[1]\n";
#
#------------------------------------------------------------------------------------
#
  my ($state,$max_conf,$file,$progname);
  my (@c,@c2,@vec,@more,@g,$count,$i,$nelem);
  ($state,$max_conf,$file,$progname)=@_;
  $nelem=$max_conf-1;
  if ($progname eq "columbus"){
    if (!-s "$file"){
      warn "wf_info: Cannot open the file: $file.\n";
    }else{
      open(FL,$file) or warn "wf_info: Cannot open the file: $file.\n";
      $count=0;
      while(<FL>){
        if (/indcsf     c/){
          $count++;
          if ($count == $state){
            $_=<FL>;
            for ($i=0;$i<=$max_conf-1;$i++){
              $_=<FL>;
              chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
              if (/csfs/){
                $nelem=@c-1;
                last;
              }else{
                @g=split(/\s+/,$_);
                $c[$i]=$g[1];
                $c2[$i]=$g[2];
                $vec[$i]=$g[$#g];
                $more[$i]=$g[3];
              }
            }
            last;
          }
        } # end if (/indcsf     c/){
      } # end while(<FL>){
      close(FL);
    } # end if (!-s $file){
  } # end if ($progname eq "columbus"){
  return \@c,\@c2,\@vec,\@more,$nelem;
}

#
#====================================================================================
#       TURBOMOLE-RELATED ROUTINES
#====================================================================================
#

sub test_turbomole_version{
#
#====================================================================================
#
# Return Turbomole version. 
#
#------------------------------------------------------------------------------------
#
  local ($tmfile,$version_def,$version,@g,$grb,$found,$ind,$iv);
  $version_def=6.6;
  $version=$version_def;
  $found=0;
  ($tmfile)=@_;
  open(INP,$tmfile) or warn " Cannot find $tmfile to read Turbomole version.\n";
  while(<INP>){
    if (/TURBOMOLE V/){
       chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
       @g=split(/\s+/,$_);
       $found++;
       last;
    }
  } 
  close(INP);
  if ($found == 0){
    $version=$version_def;
  }elsif($found != 0){
    $ind=0;
    foreach(@g){
       if (/TURBOMOLE/){
         $iv=$ind+1;
       }
       $ind++;
    }
    ($grb,$version)=split(/V/,$g[$iv]);
  }
  return $version;
}

sub read_control{
#
#====================================================================================
#
# read_control: read control (Turbomole), split in a array using $ delimiter.
# Return vector.
#
#------------------------------------------------------------------------------------
#
    local($dir,@vector,$i);
    ($dir)=@_;
    if ($dir ne ""){
      $dir=$dir."/";
    }
    open(CT,$dir."control") or warn "WARNNING: Cannot open control (Turbomole)!\n";
    $i=-1;
    while(<CT>){
      if (/\$/){
        $i++;
        $vector[$i]="$_";
      }else{
        $vector[$i]=$vector[$i].$_;
      }
      #print "$i: $vector[$i]\n";
    }
    close(CT);
    return @vector;
}

sub find_control{
#
#====================================================================================
#
# Find_control: search for key in control. Return index of the element.
# If not found, return -1.
#
#------------------------------------------------------------------------------------
#
   local (@vector,$key,@g,$i,$found);
   ($key,@vector)=@_;
   $found=-1;
   $i=0;
   foreach (@vector){
      @g=split(/\s+/,$_);
      if ($g[0]=~/\b$key\b/i){
         $found=$i;
         last;
      }
      $i++;
   }
   return $found;
}

sub get_control{
#
#====================================================================================
#
# Get_control: search for key and return value. If not found, return "not found";
#
#------------------------------------------------------------------------------------
#
   local($key,$value,@vector);
   ($key,@vector)=@_;
   $value="not found";
   $found=find_control($key,@vector);
   if ($found != -1){
     $vector[$found] =~ s/\$$key//;
     $value = $vector[$found];
     chomp($value);
     $value =~ s/^\s*//;
     $value =~ s/\s*$//;
   }
   return $value;
}

sub change_control{
#
#====================================================================================
#
# Change_control: search for key in control. Change value. If not found, add
# "key value" to the second last position.
#
#------------------------------------------------------------------------------------
#
    local (@vector,$key,$value,$found,$last);
    ($key,$value,@vector)=@_;
    $found=find_control($key,@vector);
    if ($found == -1){
       $last=$vector[-1];
       pop(@vector);
       push(@vector,"\$$key $value\n");
       push(@vector,$last);
    }else{
       $vector[$found]="\$$key $value\n";
    }
    return @vector;
}

sub delete_control{
#
#====================================================================================
#
# delete_control: delete key from vector and return vector
#
#------------------------------------------------------------------------------------
#
    local (@vector,$key,$found,$i);
    ($key,@vector)=@_;
    $found=find_control($key,@vector);
    if ($found != -1){
      delete $vector[$found];
    }
    return @vector;
}

sub print_control{
#
#====================================================================================
#
# print_control
#
#------------------------------------------------------------------------------------
#
    local (@vector);
    (@vector)=@_;
    open(CT,">control") or warn "WARNNING: Cannot write control (Turbomole)!\n";
    foreach (@vector){
      print CT "$_";
    }
    close(CT);
}

sub get_d1{
#
#====================================================================================
#
# Reads D1 diagnostic in Turbomole output file ($file).
# $method = MP2, CC2
#
#------------------------------------------------------------------------------------
#
  my ($value,$file,$method,@g);
  ($method,$file)=@_;
  $value = -1;
  open(FL,$file) or die "Cannot read $file";
  while(<FL>){
    if (/Final $method energy/){
       while(<FL>){
          if (/D1 diagnostic/){
             chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
             @g=split(/\s+/,$_);
             $value=$g[4];
             last;
          }
       }
       last;
    }
  }
  close(FL);
  return $value;
}

#
#====================================================================================
#       GAUSSIAN-RELATED ROUTINES
#====================================================================================
#

sub lookatgau{
#
#====================================================================================
#
# Look at Gaussian.
#
#------------------------------------------------------------------------------------
#
          $file = $_[0];
          open(RE,"$file") or warn "Cannot open $file!"; # I'm not sure: die or warn?
          while(<RE>){
             if (/Error termination/){
               if (/9999/){
                 # print_STDOUT("WARN: Gaussian finished in error LINK 9999.\n");
               }else{
                 print_STDOUT("ERROR: non 9999 error found in $file. \nThe message is: $_\n");
                 $status = "with ERROR";
               }
             }
          }
          close(RE);
}

#
#====================================================================================
#       GAMESS-RELATED ROUTINES
#====================================================================================
#

sub make_gamess_scr{
#
#====================================================================================
#
# Set scratch.
#
#------------------------------------------------------------------------------------
#
     my ($pathscr) = @_;
     if (! -e "$pathscr"){
        mkdir "$pathscr";
     }
 }

sub test_gamess{
#
#====================================================================================
#
# Check executable.
#
#------------------------------------------------------------------------------------
#
     my ($GAMESS, $verno, $run_gamess, $message, $scrdir, $mdle);  
     ($GAMESS, $verno, $run_gamess, $scrdir) = @_;
     $mdle = "test_gamess subroutine";
     $message = "";
     if (! -e "$GAMESS/gamess.$verno.x"){
        $message = "Cannot find GAMESS executable gamess.$verno.x in $GAMESS.\n              Check verno keyword in gamess.par.";
     }
     if (! -e "$GAMESS/$run_gamess"){
        $message = $message."\nCannot find GAMESS script $run_gamess in $GAMESS.\n              Check run_gamess keyword in gamess.par.";
     }
     if (! -e "$scrdir"){
        $message = $message."\nCannot find GAMESS scratch $scrdir.\n              Check scr keyword in gamess.par.";
     }
     if ($message ne ""){
        exit_error("$message", $mdle);
     }
 }

sub find_scr_gamess{
#
#====================================================================================
#
# Search executable.
#
#------------------------------------------------------------------------------------
#
     my ($GAMESS, $run_gamess, $value, $grb);
     ($GAMESS, $run_gamess) = @_;
     open(IN, "$GAMESS/$run_gamess") or warn "Cannot open $GAMESS/$run_gamess to read!";
     $value = "not found";
     while(<IN>){
          if (/set \s*SCR\W/){           # Look for the fisrt "set SCR"
             ($grb, $value) = split(/=/, $_);
             chomp($value);
             $value =~ s/^\s*//;
             $value =~ s/\s*$//;
             $grb = "";
             last;
          }
     }
     close(IN);
     return $value;
 }


#
#====================================================================================
#       CP2K-RELATED ROUTINES
#====================================================================================
#read_cp2kinput: read CP2K input, split in a array using \n delimiter and return vector.
 sub read_cp2kinp{
    my($dir,@vector,$ivec);

    ($dir)=@_;
    if ($dir ne ""){
      $dir=$dir."/";
    }

    open(CT,$dir."cp2k.inp") or warn "WARNNING: Cannot open CP2K input!\n";
    $ivec=-1;
    while(<CT>){
      if (/\n/){
        $ivec++;
        $vector[$ivec]="$_";
        $vector[$ivec] =~ s/^\s+|\s+$//g;
        $vector[$ivec] =~ s/\s+/ /g;
      }
    }
    close(CT);
    return @vector;
}
#---------------------------------------------------------------------------
# Get CP2K input: search for key and return value. If not found, return "not found";
sub get_cp2kinp{

   my($key,$value,@vector,$found);

   ($key,@vector)=@_;
   $value="not found";
   $found=find_cp2kinp($key,@vector);
   if ($found != -1){
     $value = $vector[$found];
   }
   return $value;
}    
#---------------------------------------------------------------------------
# find_cp2kinp: search for input line in cp2kinp. Return index of the element. If not found, return -1.
sub find_cp2kinp{
   my (@vector,$key,$i,$found);
   ($key,@vector)=@_;
   $found=-1;
   $i=0;
      foreach (@vector){
        if ($vector[$i]=~/$key/i){
          $found=$i;
          last;
          }
       $i++;
       }
   return $found;
}
#---------------------------------------------------------------------------
# find_cp2kinp: search for input line in cp2kinp. Return index of the element. If not found, return -1.
sub find_cp2kinp_full{
   local (@vector,$key,@g,$i,$found);
   ($key,@vector)=@_;
   $found=-1;
   $i=0;
      foreach (@vector){
        if ($vector[$i]=~/^$key$/i){
          $found=$i;
          last;
          }
       $i++;
       }
   return $found;
}
#---------------------------------------------------------------------------
# change_cp2kinput: search for key in CP2K input. Find value and change it if necessary. If not found, exit.
sub change_cp2kinp{

    my (@vector,$key,$value,$found,$last,$keyforfind,$string,$foundkey,$oldinput);

    ($key,$value,@vector,$mdle)=@_;
    $keyforfind="$key"." "."$value";
    $found=find_cp2kinp_full($keyforfind,@vector);
    if ($found == -1){
            print "WARNING: $key is NOT set to $value in CP2K input, but requested.\n";
            $string="$key"." ";
            $foundkey=find_cp2kinp($string,@vector);
            if ($foundkey != -1){
               $string="$key"." "."$value";
               $oldinput=$vector[$foundkey];
               print "WARNING: Thus resetting OLD INPUT $vector[$foundkey] to NEW INPUT $string.\n";
               $vector[$foundkey]="$string";
               system("sed -i '/$oldinput/c\ $string' cp2k.inp");
            }else{
               exit_error("$key keyword / section missing in CP2K input.",$mdle);
            }
    }
    return @vector;
}
#------------------------------------------------------------------------------
# rm_cp2kinput: rm CP2K input section
sub rm_cp2kinp{

    my(@vector,$key,$value,$found,$ii,$sizevector,$ends_to_count,$starts_found,$sizevector,$length);

    ($key,@vector)=@_;
    $found=find_cp2kinp($key,@vector);
    if ($found == -1){
        print "WARNING: Cannot remove $key from CP2K input, because $key not found.\n";
    }else{
      $index=$found;
      $ends_to_count=0;
      $starts_found=0;
      $sizevector=@vector;
      for ($ii=$index;$ii<=$sizevector;$ii++){
       $first=substr($vector[$ii],0,1);
       if ($first eq "\&"){
         $starts_found=$starts_found+1;
         $value=substr($vector[$ii],0,4);
         if ($value eq "\&END"){
          $starts_found=$starts_found-1; #&END should not be counted as new keyword section
          $ends_to_count=$ends_to_count+1;
         if ($ends_to_count == $starts_found){
          $secondindex=$ii;
          last;
         }
         }
       }
       }
      $length=1+$secondindex-$index;
      splice(@vector,$index,$length);

      $index=$index+1;
      $secondindex=$secondindex+1;
      system("sed -i '$index,$secondindex d' cp2k.inp"); # only file is changed, not array!
    }
    return @vector;
}
#---------------------------------------------------------------------
# add_cp2kinput: add CP2K input section
sub add_cp2kinp{

    my(@vector,@tmparray,$key,$sizevector,$ii,$jj,$kk,$string,$sizeend,$stringstart,$endstring);

    ($key,@vector)=@_;
    @tmparray=split(/\//,$key);
    $sizevector=@tmparray-1;

    $string="\&$tmparray[1]";
    $sizeend=$sizevector-1;
    $endstring="\&END $tmparray[$sizeend]";

    for ($ii=1;$ii<=$sizevector;$ii++){
         $jj=$ii+1;
         $kk=$sizevector-$ii-1;
     if ($jj<$sizevector){
        $string="$string"."\\n \&$tmparray[$jj]";
        if ($kk>0){
           $endstring="$endstring"."\\n \&END $tmparray[$kk]";
        }
     }else{
     $string="$string"."\\n $tmparray[$jj]";
    }
    }

    $stringstart="\&$tmparray[0]";
    $endstring="$string"."$endstring";
    system("sed -i '/$stringstart/a\ $endstring' cp2k.inp"); #only file is changed, not array!
}
#---------------------------------------------------------------------
# add_cp2kinput: add CP2K keyword
sub add_cp2kkey{

    my(@vector,@tmparray,$key,$string,$stringstart);

    ($key,@vector)=@_;
    @tmparray=split(/\//,$key);

    $stringstart="\&$tmparray[0]";
    $string="$tmparray[1]";

    system("sed -i '/$stringstart/a\ $string' cp2k.inp"); #only file is changed, not array!
}
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#====================================================================================
#                               END SUBROUTINES
#====================================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
