# Fromage dynamic test cases

    Test  initial_state	   singlet	 triplet   single_state   nactype   SH method    coupling	                          description
    -------------------------------------------------------------------------------------------------------------------------------------------------------------------
     1	     2	             2	           0	        0	    ktdc      FSSH       []                                     test default values  
     2	     2	             3	           0	      	1	    ktdc      FSSH       [1 2 2 3]                              test single state grad with approx NAC for FSSH 
     3	     2	             3	           0	      	1	    nac       FSSH       [1 2 2 3]                              test single state grad with exact NAC for FSSH
     4	     2	             3	           0	      	0	    n/a	      GSH        n/a                                    test all state grad for GSH
     5	     3	             2	           2	      	1	    nac       FSSH       [1 2 3 4 1 3 2 3 1 4 2 4]              test nac and soc calculation
     6	     2	             3	           2	      	1	    ktdc      FSSH       [1 2 2 3 4 5 1 4 1 5 2 4 2 5 3 4 3 5]  test FSSH with SOC
     7	     2	             3	           2	      	0	    n/a	      GSH        n/a                                    test GSH with SOC
