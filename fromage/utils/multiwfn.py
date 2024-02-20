
def calc_mwfn_charges(cmd = "Multiwfn_noGUI", mwfn_in = "mwfn.in", molden_in = "molden.input", charge_type="Mulliken", nprocs =None ):
    """post-process xtb wavefunction (.molden format) to type of charge"""
    
    mwfn_inputs = {
    "mulliken": "7\n5\n1\ny\n0\n0\n\nq",
    "lowdin": "7\n6\n1\ny\n0\n0\n\nq",
    "hirshfeld": "7\n11\n1\ny\n0\n0\n\nq",
    "resp": "7\n18\n1\ny\n0\n0\n\nq",
    "chelpg": "7\n12\n1\ny\n0\n0\n\nq",
    "cm5": "7\n16\n1\ny\n0\n\nq",
    }
    #'molden_file = molden_file
    exc_inp = mwfn_inputs[charge_type.lower()]
    
    with open(mwfn_in, "w") as mwfn_in:
        mwfn_in.write(exc_inp)
    
    if not nprocs == None:
        # clean up previous charges
        if os.path.exists('molden.chg'):
            subprocess.call("rm molden.chg", shell=True)

        # this needs modifying for safety against injection with shell=True 
        result = subprocess.run("{} {} -nt {} < mwfn.in".format(cmd, molden_in, nprocs), shell=True)
       
        print("Multiwfn ran successfully")
    else:
        raise ValueError("Mutlwifn requires specification of number of processers")        
   
    
    mwfn_charges = read_mwfn_charges("molden.chg")

    return mwfn_charges