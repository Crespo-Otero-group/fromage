%chk=naphthalene222.clust.chk
%mem=2GB
%nproc=4
#wb97xd 6-311++g** opt td(nstates=1,root=1)

naphthalene222.clust

0 1
     C  -1.032790  -3.082595  -9.226208
     C  -1.410845  -3.953473 -10.244068
     H  -0.221479  -3.368330  -8.539660
     C  -1.694137  -1.828081  -9.066557
     C  -2.463407  -3.602963 -11.143879
     H  -0.911067  -4.924554 -10.374742
     C  -2.722477  -1.454857  -9.927332
     H  -1.397175  -1.146999  -8.253729
     C  -2.873226  -4.480042 -12.194287
     C  -3.132296  -2.331935 -10.977741
     H  -3.238263  -0.489802  -9.810230
     C  -3.901566  -4.106818 -13.055063
     H  -2.357440  -5.445097 -12.311390
     C  -4.184858  -1.981425 -11.877552
     H  -4.198527  -4.787899 -13.867891
     C  -4.562913  -2.852303 -12.895412
     H  -4.684636  -1.010344 -11.746878
     H  -5.374224  -2.566568 -13.581959
     C  -2.708435  -5.819738  -5.539255
     C  -3.086482  -4.948868  -6.557129
     C  -3.369749  -7.074278  -5.379686
     H  -1.897122  -5.534025  -4.852700
     C  -4.139063  -5.299377  -7.456915
     H  -2.586709  -3.977788  -6.687818
     C  -4.398161  -7.447486  -6.240373
     H  -3.072812  -7.755341  -4.566828
     C  -4.808003  -6.570380  -7.290739
     C  -4.548853  -4.422335  -8.507373
     H  -4.913926  -8.412563  -6.123348
     C  -5.860477  -6.920920  -8.190667
     C  -5.577218  -4.795524  -9.368149
     H  -4.033082  -3.457260  -8.624421
     C  -6.238516  -6.050074  -9.208572
     H  -6.360281  -7.891980  -8.059919
     H  -5.874174  -4.114445 -10.180972
     H  -7.049895  -6.335783  -9.895052
     C   0.292188  -7.915428  -4.506090
     C  -0.082585  -8.788075  -5.522539
     H  -0.209624  -6.945947  -4.372777
     C   1.341237  -8.266944  -3.604910
     C   0.578086 -10.042703  -5.680517
     H  -0.895701  -8.503642  -6.207594
     C   2.010063  -9.537753  -3.768973
     C   1.744705  -7.389667  -2.553405
     C   1.606595 -10.415030  -4.820478
     H   0.278928 -10.723262  -6.492899
     C   3.059111  -9.889269  -2.867793
     C   2.773214  -7.761995  -1.693365
     H   1.232992  -6.422936  -2.437050
     H   2.118308 -11.381762  -4.936834
     C   3.433885  -9.016623  -1.851343
     H   3.560924 -10.858751  -3.001106
     H   3.072372  -7.081436  -0.880984
     H   4.247001  -9.301056  -1.166289
     C  -5.814654  -3.602938  -3.770049
     C  -6.224468  -4.480007  -4.820464
     C  -4.762163  -3.953477  -2.870144
     C  -6.483647  -2.331961  -3.603834
     C  -7.252904  -4.106835  -5.681151
     H  -5.708718  -5.445105  -4.937435
     C  -4.384130  -3.082639  -1.852227
     H  -4.262364  -4.924536  -3.000909
     C  -7.536137  -1.981421  -4.503739
     C  -6.073832  -1.454891  -2.553419
     H  -7.549836  -4.787895  -6.494004
     C  -7.914171  -2.852259  -5.521656
     H  -3.572750  -3.368326  -1.165739
     C  -5.045396  -1.828063  -1.692732
     H  -8.035936  -1.010362  -4.372973
     H  -6.589583  -0.489793  -2.436448
     H  -8.725551  -2.566572  -6.208144
     H  -4.748464  -1.147003  -0.879879
     C   1.609257  -4.482799  -4.819422
     C   0.581235  -4.111585  -5.679723
     H   2.120061  -5.449596  -4.936849
     C   2.010649  -3.603397  -3.768714
     H   0.278575  -4.789182  -6.493873
     C  -0.073239  -2.855476  -5.511902
     C   1.344059  -2.331749  -3.601645
     C   3.059515  -3.954741  -2.867620
     C   0.295187  -1.976765  -4.500120
     H  -0.885388  -2.572259  -6.198821
     C   1.745544  -1.453114  -2.551350
     C   3.434156  -3.082018  -1.851256
     H   3.561068  -4.924210  -3.001198
     H  -0.210802  -1.008668  -4.372419
     C   2.774030  -1.827215  -1.692561
     H   1.232274  -0.486950  -2.437066
     H   4.247197  -3.366184  -1.166095
     H   3.072383  -1.146127  -0.880681
     C  -1.758283  -6.049141  -1.835508
     C  -1.383484  -6.921824  -0.819111
     C  -1.097557  -4.794554  -1.993567
     H  -2.571329  -6.333603  -2.520633
     C  -0.334354  -6.570333   0.081958
     H  -1.885277  -7.891322  -0.685867
     C  -0.069069  -4.422190  -1.133530
     H  -1.396724  -4.113991  -2.805950
     C   0.334430  -5.299505  -0.082076
     C   0.069095  -7.447587   1.133501
     H   0.442653  -3.455476  -1.249931
     C   1.383454  -4.947986   0.819134
     C   1.097522  -7.075237   1.993640
     H  -0.442641  -8.414295   1.249936
     C   1.758238  -5.820644   1.835565
     H   1.885270  -3.978502   0.685839
     H   1.396710  -7.755820   2.805984
     H   2.571346  -5.536185   2.520618
     C  -3.433935  -8.788126   1.851451
     C  -3.059126  -7.915431   2.867831
     C  -2.773164 -10.042679   1.693298
     H  -4.246973  -8.503638   1.166328
     C  -2.010019  -8.266924   3.768929
     H  -3.560922  -6.945930   3.001056
     C  -1.744758 -10.415064   2.553432
     H  -3.072361 -10.723265   0.880951
     C  -1.341279  -9.537773   3.604925
     C  -1.606540  -7.389634   4.820422
     H  -1.233013 -11.381755   2.436952
     C  -0.292172  -9.889267   4.506024
     C  -0.578134  -7.762019   5.680557
     H  -2.118285  -6.422943   4.936903
     C   0.082637  -9.016572   5.522403
     H   0.209624 -10.858767   4.372798
     H  -0.278937  -7.081433   6.492903
     H   0.895675  -9.301060   6.207526
     C  -5.109536  -6.049154   5.538289
     C  -4.734752  -6.921812   6.554721
     C  -4.448819  -4.794561   5.380215
     H  -5.922644  -6.333613   4.853236
     C  -3.685727  -6.570293   7.455930
     H  -5.236567  -7.891296   6.688016
     C  -3.420392  -4.422211   6.240353
     H  -4.748008  -4.113978   4.567870
     C  -3.016944  -5.299464   7.291896
     C  -3.282230  -7.447606   8.507386
     H  -2.908656  -3.455503   6.123918
     C  -1.967813  -4.947974   8.192965
     C  -2.253741  -7.075244   9.367421
     H  -3.793950  -8.414322   8.623786
     C  -1.593015  -5.820657   9.209363
     H  -1.466020  -3.978475   8.059722
     H  -1.954574  -7.755807  10.179804
     H  -0.779969  -5.536195   9.894487
     C  -7.735387  -3.082595   5.521529
     C  -8.113442  -3.953473   4.503669
     H  -6.924076  -3.368330   6.208077
     C  -8.396734  -1.828081   5.681180
     C  -9.166004  -3.602963   3.603858
     H  -7.613664  -4.924554   4.372995
     C  -9.425074  -1.454857   4.820405
     H  -8.099773  -1.146999   6.494008
     C  -9.575824  -4.480042   2.553450
     C  -9.834894  -2.331935   3.769996
     H  -9.940860  -0.489802   4.937507
     C -10.604163  -4.106818   1.692674
     H  -9.060038  -5.445097   2.436348
     C -10.887456  -1.981425   2.870185
     H -10.901125  -4.787899   0.879846
     C -11.265510  -2.852303   1.852325
     H -11.387233  -1.010344   3.000859
     H -12.076822  -2.566568   1.165778
     C  -6.059784  -5.819724   1.834689
     C  -6.437824  -4.948878   0.816785
     C  -6.721083  -7.074274   1.994266
     H  -5.248406  -5.534015   2.521169
     C  -7.490296  -5.299418  -0.083144
     H  -5.938019  -3.977818   0.686036
     C  -7.749447  -7.447463   1.133490
     H  -6.424126  -7.755353   2.807089
     C  -8.159237  -6.570420   0.083033
     C  -7.900140  -4.422312  -1.133509
     H  -8.265218  -8.412538   1.250538
     C  -9.211817  -6.920929  -0.816754
     C  -8.928551  -4.795520  -1.994196
     H  -7.384374  -3.457235  -1.250534
     C  -9.589866  -6.050060  -1.834628
     H  -9.711590  -7.892009  -0.686064
     H  -9.225488  -4.114457  -2.807055
     H -10.401179  -6.335773  -2.521183
     C  -1.744758  -4.480164   2.553432
     C  -2.773163  -4.107779   1.693297
     H  -1.233013  -5.446855   2.436951
     C  -1.341279  -3.602874   3.604925
     H  -3.072361  -4.788365   0.880951
     C  -3.433935  -2.853226   1.851451
     C  -2.010019  -2.332025   3.768929
     C  -0.292172  -3.954367   4.506024
     C  -3.059126  -1.980531   2.867831
     H  -4.246973  -2.568738   1.166328
     C  -1.606540  -1.454735   4.820422
     C   0.082637  -3.081672   5.522403
     H   0.209624  -4.923868   4.372798
     H  -3.560922  -1.011030   3.001056
     C  -0.578134  -1.827119   5.680557
     H  -2.118285  -0.488043   4.936903
     H   0.895675  -3.366160   6.207526
     H  -0.278937  -1.146533   6.492903
     C  -6.073832  -7.389791  -2.553419
     C  -6.483647  -8.266860  -3.603834
     C  -5.045396  -7.762963  -1.692732
     H  -6.589583  -6.424693  -2.436449
     C  -7.536137  -7.916321  -4.503739
     C  -5.814654  -9.537837  -3.770049
     C  -4.384130  -9.017539  -1.852227
     H  -4.748464  -7.081903  -0.879879
     C  -7.914171  -8.787159  -5.521656
     H  -8.035936  -6.945262  -4.372974
     C  -6.224468 -10.414907  -4.820463
     C  -4.762163  -9.888377  -2.870144
     H  -3.572750  -9.303226  -1.165739
     C  -7.252904 -10.041735  -5.681151
     H  -8.725551  -8.501472  -6.208144
     H  -5.708718 -11.380005  -4.937434
     H  -4.262364 -10.859436  -3.000909
     H  -7.549836 -10.722795  -6.494004
     C   1.593061  -0.114255  -9.209448
     C   1.967845  -0.986913  -8.193016
     C   2.253778   1.140339  -9.367522
     H   0.779954  -0.398714  -9.894501
     C   3.016870  -0.635394  -7.291808
     H   1.466030  -1.956396  -8.059721
     C   3.282205   1.512688  -8.507384
     H   1.954590   1.820921 -10.179867
     C   3.685654   0.635435  -7.455841
     C   3.420369  -1.512707  -6.240352
     H   3.793941   2.479397  -8.623819
     C   4.734784   0.986926  -6.554772
     C   4.448857  -1.140345  -5.380316
     H   2.908647  -2.479422  -6.123951
     C   5.109583   0.114243  -5.538374
     H   5.236577   1.956425  -6.688015
     H   4.748024  -1.820907  -4.567933
     H   5.922629   0.398705  -4.853250
     C  -0.082585   3.081725  -5.522539
     C   0.292188   3.954371  -4.506090
     C   0.578086   1.827097  -5.680517
     H  -0.895701   3.366158  -6.207594
     H  -0.209624   4.923853  -4.372777
     C   1.341237   3.602855  -3.604910
     C   1.606595   1.454769  -4.820478
     H   0.278928   1.146538  -6.492899
     C   2.010063   2.332046  -3.768973
     C   1.744705   4.480132  -2.553405
     H   2.118308   0.488038  -4.936834
     C   3.059111   1.980530  -2.867793
     C   2.773214   4.107805  -1.693365
     H   1.232992   5.446864  -2.437050
     C   3.433885   2.853177  -1.851343
     H   3.560924   1.011049  -3.001106
     H   3.072372   4.788363  -0.880984
     H   4.247001   2.568744  -1.166289
     C  -2.708435   0.115162  -5.539255
     C  -3.086483   0.986031  -6.557129
     C  -3.369749  -1.139379  -5.379687
     H  -1.897122   0.400875  -4.852700
     C  -4.139063   0.635522  -7.456915
     H  -2.586710   1.957111  -6.687819
     C  -4.398161  -1.512586  -6.240373
     H  -3.072812  -1.820441  -4.566828
     C  -4.808004  -0.635480  -7.290739
     C  -4.548853   1.512564  -8.507373
     H  -4.913926  -2.477663  -6.123349
     C  -5.860477  -0.986020  -8.190667
     C  -5.577218   1.139376  -9.368149
     H  -4.033082   2.477640  -8.624421
     C  -6.238516  -0.115174  -9.208572
     H  -6.360281  -1.957081  -8.059919
     H  -5.874174   1.820455 -10.180972
     H  -7.049895  -0.400883  -9.895052
     C  -1.032790   2.852305  -9.226208
     C  -1.410845   1.981427 -10.244068
     C  -1.694137   4.106819  -9.066557
     H  -0.221479   2.566570  -8.539661
     C  -2.463407   2.331937 -11.143879
     H  -0.911067   1.010346 -10.374742
     C  -2.722477   4.480043  -9.927332
     H  -1.397175   4.787901  -8.253729
     C  -3.132296   3.602965 -10.977741
     C  -2.873226   1.454858 -12.194287
     H  -3.238263   5.445098  -9.810230
     C  -4.184858   3.953474 -11.877552
     C  -3.901566   1.828082 -13.055063
     H  -2.357440   0.489803 -12.311389
     C  -4.562913   3.082597 -12.895412
     H  -4.684636   4.924555 -11.746878
     H  -4.198527   1.147001 -13.867891
     H  -5.374224   3.368332 -13.581959
     C  -2.708435   6.050061  -5.539255
     C  -3.086482   6.920931  -6.557129
     C  -3.369749   4.795521  -5.379686
     H  -1.897122   6.335774  -4.852700
     C  -4.139063   6.570423  -7.456915
     H  -2.586709   7.892011  -6.687818
     C  -4.398161   4.422314  -6.240373
     H  -3.072812   4.114458  -4.566828
     C  -4.808003   5.299420  -7.290739
     C  -4.548853   7.447464  -8.507373
     H  -4.913926   3.457236  -6.123348
     C  -5.860477   4.948880  -8.190667
     C  -5.577218   7.074276  -9.368149
     H  -4.033082   8.412540  -8.624421
     C  -6.238516   5.819726  -9.208572
     H  -6.360281   3.977819  -8.059919
     H  -5.874174   7.755355 -10.180972
     H  -7.049895   5.534017  -9.895052
     C   1.967845   4.947987  -8.193016
     C   1.593061   5.820645  -9.209448
     C   3.016870   5.299507  -7.291807
     H   1.466030   3.978504  -8.059721
     C   2.253778   7.075239  -9.367522
     H   0.779954   5.536186  -9.894501
     C   3.685654   6.570335  -7.455841
     C   3.420367   4.422193  -6.240351
     C   3.282205   7.447588  -8.507384
     H   1.954590   7.755821 -10.179867
     C   4.734784   6.921826  -6.554772
     C   4.448856   4.794556  -5.380316
     H   2.908647   3.455478  -6.123951
     H   3.793941   8.414297  -8.623819
     C   5.109583   6.049143  -5.538374
     H   5.236577   7.891324  -6.688015
     H   4.748024   4.113992  -4.567933
     H   5.922628   6.333605  -4.853250
     C  -7.536137   3.953478  -4.503739
     C  -7.914171   3.082641  -5.521656
     H  -8.035936   4.924538  -4.372974
     C  -6.483647   3.602939  -3.603834
     C  -7.252904   1.828065  -5.681151
     H  -8.725551   3.368328  -6.208144
     C  -5.814654   2.331962  -3.770049
     C  -6.073832   4.480009  -2.553419
     C  -6.224468   1.454893  -4.820463
     H  -7.549836   1.147004  -6.494004
     C  -4.762163   1.981423  -2.870144
     C  -5.045396   4.106837  -1.692732
     H  -6.589583   5.445107  -2.436449
     H  -5.708718   0.489795  -4.937434
     C  -4.384130   2.852261  -1.852227
     H  -4.262364   1.010363  -3.000909
     H  -4.748464   4.787897  -0.879879
     H  -3.572750   2.566574  -1.165739
     C   1.609257   7.387001  -4.819422
     C   0.581235   7.758215  -5.679723
     H   2.120061   6.420204  -4.936849
     C   2.010649   8.266402  -3.768714
     H   0.278575   7.080617  -6.493873
     C  -0.073239   9.014323  -5.511902
     C   1.344059   9.538051  -3.601645
     C   3.059515   7.915059  -2.867620
     C   0.295187   9.893034  -4.500120
     H  -0.885388   9.297541  -6.198821
     C   1.745544  10.416686  -2.551350
     C   3.434156   8.787781  -1.851256
     H   3.561068   6.945590  -3.001198
     H  -0.210802  10.861132  -4.372419
     C   2.774030  10.042585  -1.692561
     H   1.232274  11.382850  -2.437066
     H   4.247197   8.503615  -1.166095
     H   3.072383  10.723672  -0.880681
     H  -5.708718   6.424695  -4.937435
     C  -6.224468   7.389793  -4.820464
     C  -5.814654   8.266862  -3.770049
     C  -7.252904   7.762964  -5.681151
     C  -4.762163   7.916323  -2.870144
     C  -6.483647   9.537839  -3.603834
     H  -7.549836   7.081904  -6.494004
     C  -7.914171   9.017541  -5.521656
     C  -4.384130   8.787160  -1.852227
     H  -4.262364   6.945263  -3.000909
     C  -7.536137   9.888378  -4.503739
     C  -6.073832  10.414909  -2.553419
     H  -8.725551   9.303227  -6.208144
     H  -3.572750   8.501474  -1.165739
     C  -5.045396  10.041737  -1.692732
     H  -8.035936  10.859438  -4.372973
     H  -6.589583  11.380006  -2.436448
     H  -4.748464  10.722797  -0.879879
     C  -5.109536  -0.114255   5.538289
     C  -4.734752  -0.986913   6.554721
     C  -4.448819   1.140339   5.380215
     H  -5.922644  -0.398714   4.853236
     C  -3.685727  -0.635394   7.455929
     H  -5.236568  -1.956396   6.688016
     C  -3.420392   1.512688   6.240353
     H  -4.748008   1.820921   4.567870
     C  -3.016944   0.635435   7.291896
     C  -3.282229  -1.512707   8.507385
     H  -2.908656   2.479397   6.123918
     C  -1.967813   0.986926   8.192965
     C  -2.253741  -1.140345   9.367421
     H  -3.793950  -2.479422   8.623786
     C  -1.593015   0.114243   9.209363
     H  -1.466020   1.956425   8.059722
     H  -1.954574  -1.820907  10.179804
     H  -0.779969   0.398705   9.894487
     C  -1.758283   5.820659  -1.835508
     C  -1.383484   4.947975  -0.819111
     C  -1.097557   7.075246  -1.993567
     H  -2.571329   5.536196  -2.520633
     C  -0.334354   5.299467   0.081958
     H  -1.885277   3.978477  -0.685867
     C  -0.069069   7.447609  -1.133530
     H  -1.396724   7.755809  -2.805950
     C   0.334430   6.570295  -0.082076
     C   0.069095   4.422213   1.133501
     H   0.442653   8.414324  -1.249931
     C   1.383454   6.921814   0.819134
     C   1.097522   4.794563   1.993640
     H  -0.442641   3.455505   1.249936
     C   1.758238   6.049156   1.835565
     H   1.885270   7.891297   0.685839
     H   1.396710   4.113980   2.805984
     H   2.571346   6.333615   2.520618
     C  -3.433935   3.081674   1.851451
     C  -3.059126   3.954369   2.867831
     C  -2.773164   1.827121   1.693298
     H  -4.246973   3.366162   1.166328
     C  -2.010019   3.602875   3.768929
     H  -3.560922   4.923869   3.001056
     C  -1.744758   1.454736   2.553432
     H  -3.072361   1.146535   0.880951
     C  -1.341279   2.332026   3.604925
     C  -1.606540   4.480165   4.820422
     H  -1.233013   0.488045   2.436952
     C  -0.292172   1.980533   4.506024
     C  -0.578134   4.107781   5.680557
     H  -2.118285   5.446856   4.936903
     C   0.082637   2.853227   5.522403
     H   0.209624   1.011032   4.372798
     H  -0.278937   4.788367   6.492903
     H   0.895675   2.568739   6.207526
     C  -5.109536   5.820645   5.538289
     C  -4.734752   4.947987   6.554721
     C  -4.448819   7.075239   5.380215
     H  -5.922644   5.536186   4.853236
     C  -3.685727   5.299507   7.455930
     H  -5.236567   3.978504   6.688016
     C  -3.420392   7.447588   6.240353
     H  -4.748008   7.755821   4.567870
     C  -3.016944   6.570335   7.291896
     C  -3.282230   4.422193   8.507386
     H  -2.908656   8.414297   6.123918
     C  -1.967813   6.921826   8.192965
     C  -2.253741   4.794556   9.367421
     H  -3.793950   3.455478   8.623786
     C  -1.593015   6.049143   9.209363
     H  -1.466020   7.891324   8.059722
     H  -1.954574   4.113992  10.179804
     H  -0.779969   6.333605   9.894487
     C  -6.059784   0.115176   1.834689
     C  -6.437824   0.986022   0.816785
     C  -6.721082  -1.139375   1.994266
     H  -5.248406   0.400885   2.521169
     C  -7.490296   0.635481  -0.083144
     H  -5.938019   1.957082   0.686036
     C  -7.749447  -1.512563   1.133490
     H  -6.424126  -1.820453   2.807089
     C  -8.159238  -0.635521   0.083033
     C  -7.900140   1.512588  -1.133509
     H  -8.265218  -2.477638   1.250538
     C  -9.211817  -0.986030  -0.816754
     C  -8.928551   1.139380  -1.994196
     H  -7.384374   2.477665  -1.250534
     C  -9.589866  -0.115160  -1.834627
     H  -9.711590  -1.957110  -0.686064
     H  -9.225488   1.820443  -2.807055
     H -10.401179  -0.400873  -2.521182
     C  -6.059784   6.050076   1.834689
     C  -6.437824   6.920922   0.816785
     C  -6.721083   4.795525   1.994266
     H  -5.248406   6.335785   2.521169
     C  -7.490296   6.570381  -0.083144
     H  -5.938019   7.891982   0.686036
     C  -7.749447   4.422337   1.133490
     H  -6.424126   4.114447   2.807089
     C  -8.159237   5.299379   0.083033
     C  -7.900140   7.447488  -1.133509
     H  -8.265218   3.457261   1.250538
     C  -9.211817   4.948871  -0.816754
     C  -8.928551   7.074280  -1.994196
     H  -7.384374   8.412565  -1.250534
     C  -9.589866   5.819740  -1.834628
     H  -9.711590   3.977790  -0.686064
     H  -9.225488   7.755343  -2.807055
     H -10.401179   5.534027  -2.521183
     C  -7.735387   2.852305   5.521529
     C  -8.113442   1.981427   4.503669
     C  -8.396734   4.106819   5.681180
     H  -6.924076   2.566570   6.208076
     C  -9.166004   2.331937   3.603858
     H  -7.613664   1.010346   4.372995
     C  -9.425074   4.480043   4.820405
     H  -8.099773   4.787901   6.494008
     C  -9.834893   3.602965   3.769996
     C  -9.575824   1.454858   2.553450
     H  -9.940860   5.445098   4.937507
     C -10.887456   3.953474   2.870185
     C -10.604164   1.828082   1.692675
     H  -9.060038   0.489803   2.436348
     C -11.265510   3.082597   1.852325
     H -11.387233   4.924555   3.000859
     H -10.901125   1.147001   0.879846
     H -12.076822   3.368332   1.165778
     C  -1.744758   7.389636   2.553432
     C  -2.773163   7.762021   1.693297
     H  -1.233013   6.422945   2.436951
     C  -1.341279   8.266926   3.604925
     H  -3.072361   7.081434   0.880951
     C  -3.433935   9.016574   1.851451
     C  -2.010019   9.537775   3.768929
     C  -0.292172   7.915432   4.506024
     C  -3.059126   9.889269   2.867831
     H  -4.246973   9.301062   1.166328
     C  -1.606540  10.415065   4.820422
     C   0.082637   8.788127   5.522403
     H   0.209624   6.945932   4.372798
     H  -3.560922  10.858769   3.001056
     C  -0.578134  10.042680   5.680557
     H  -2.118285  11.381756   4.936903
     H   0.895675   8.503639   6.207526
     H  -0.278937  10.723267   6.492903
     C   5.109583  -5.820657  -5.538374
     C   4.734784  -4.947974  -6.554772
     C   4.448856  -7.075244  -5.380316
     H   5.922628  -5.536195  -4.853250
     C   3.685654  -5.299464  -7.455841
     H   5.236577  -3.978475  -6.688015
     C   3.420367  -7.447606  -6.240351
     H   4.748024  -7.755807  -4.567933
     C   3.016870  -6.570293  -7.291807
     C   3.282205  -4.422211  -8.507384
     H   2.908647  -8.414322  -6.123951
     C   1.967845  -6.921812  -8.193016
     C   2.253778  -4.794561  -9.367522
     H   3.793941  -3.455503  -8.623819
     C   1.593061  -6.049154  -9.209448
     H   1.466030  -7.891296  -8.059721
     H   1.954590  -4.113978 -10.179867
     H   0.779954  -6.333613  -9.894501
     C   6.059734  -6.050060  -1.834628
     C   6.437783  -6.920929  -0.816754
     C   6.721049  -4.795520  -1.994196
     H   5.248421  -6.335773  -2.521183
     C   7.490363  -6.570420   0.083033
     H   5.938010  -7.892009  -0.686064
     C   7.749460  -4.422312  -1.133509
     H   6.424112  -4.114457  -2.807055
     C   8.159304  -5.299418  -0.083144
     C   7.900153  -7.447463   1.133490
     H   8.265226  -3.457235  -1.250534
     C   9.211776  -4.948878   0.816785
     C   8.928517  -7.074274   1.994266
     H   7.384382  -8.412538   1.250538
     C   9.589816  -5.819724   1.834689
     H   9.711581  -3.977818   0.686036
     H   9.225474  -7.755353   2.807089
     H  10.401194  -5.534015   2.521169
     C   4.384090  -8.787203   1.852325
     C   4.762144  -7.916325   2.870185
     C   5.045437 -10.041717   1.692675
     H   3.572778  -8.501468   1.165778
     C   5.814707  -8.266835   3.769996
     H   4.262367  -6.945244   3.000859
     C   6.073776 -10.414942   2.553450
     H   4.748475 -10.722799   0.879846
     C   6.483596  -9.537863   3.603858
     C   6.224526  -7.389756   4.820405
     H   6.589562 -11.379996   2.436348
     C   7.536158  -9.888373   4.503669
     C   7.252866  -7.762980   5.681180
     H   5.708740  -6.424702   4.937507
     C   7.914213  -9.017495   5.521529
     H   8.035936 -10.859454   4.372995
     H   7.549827  -7.081899   6.494008
     H   8.725524  -9.303230   6.208076
     C   2.708486  -6.050074   5.539165
     C   3.086526  -6.920920   6.557070
     C   3.369785  -4.795524   5.379588
     H   1.897108  -6.335783   4.852685
     C   4.138999  -6.570380   7.456998
     H   2.586722  -7.891980   6.687818
     C   4.398150  -4.422335   6.240364
     H   3.072828  -4.114445   4.566765
     C   4.807940  -5.299377   7.290822
     C   4.548842  -7.447486   8.507364
     H   4.913921  -3.457260   6.123316
     C   5.860520  -4.948868   8.190608
     C   5.577253  -7.074278   9.368051
     H   4.033077  -8.412563   8.624389
     C   6.238568  -5.819738   9.208482
     H   6.360293  -3.977788   8.059919
     H   5.874190  -7.755341  10.180909
     H   7.049881  -5.534025   9.895037
     C   7.536158  -3.953473   4.503669
     C   7.914213  -3.082595   5.521529
     C   6.483596  -3.602963   3.603858
     H   8.035936  -4.924554   4.372995
     H   8.725524  -3.368330   6.208077
     C   7.252866  -1.828081   5.681180
     C   6.073776  -4.480042   2.553450
     C   5.814706  -2.331935   3.769996
     C   6.224526  -1.454857   4.820405
     H   7.549827  -1.146999   6.494008
     C   5.045437  -4.106818   1.692674
     H   6.589563  -5.445097   2.436348
     C   4.762144  -1.981425   2.870185
     H   5.708740  -0.489802   4.937507
     H   4.748475  -4.787899   0.879846
     C   4.384090  -2.852303   1.852325
     H   4.262367  -1.010344   3.000859
     H   3.572778  -2.566568   1.165778
     H   0.442653  -9.390376  -1.249932
     C  -0.069069 -10.357091  -1.133531
     C   0.334430 -11.234405  -0.082076
     C  -1.097557 -10.729454  -1.993567
     C   1.383454 -10.882886   0.819134
     C  -0.334354 -12.505233   0.081958
     C  -1.758283 -11.984041  -1.835508
     H  -1.396724 -10.048891  -2.805950
     C   1.758238 -11.755544   1.835565
     H   1.885270  -9.913402   0.685839
     C  -1.383484 -12.856724  -0.819111
     C   0.069095 -13.382487   1.133501
     H  -2.571329 -12.268503  -2.520633
     C   1.097522 -13.010137   1.993640
     H   2.571346 -11.471084   2.520618
     H  -1.885277 -13.826223  -0.685867
     H  -0.442641 -14.349195   1.249936
     H   1.396710 -13.690719   2.805984
     H   1.397166  -4.787895   8.253733
     C   1.694099  -4.106835   9.066586
     C   2.722535  -4.480007   9.927273
     C   1.032832  -2.852259   9.226081
     C   3.132349  -3.602938  10.977688
     H   3.238285  -5.445105   9.810302
     C   1.410866  -1.981421  10.243998
     H   0.221452  -2.566572   8.539593
     C   4.184839  -3.953477  11.877593
     C   2.463356  -2.331961  11.143903
     H   0.911066  -1.010362  10.374764
     C   4.562873  -3.082639  12.895510
     H   4.684639  -4.924536  11.746828
     C   2.873170  -1.454891  12.194318
     H   5.374253  -3.368326  13.581998
     C   3.901606  -1.828063  13.055005
     H   2.357420  -0.489793  12.311289
     H   4.198539  -1.147003  13.867858
     C   7.735429  -2.852259  -5.521656
     C   8.396696  -4.106835  -5.681151
     C   8.113463  -1.981421  -4.503739
     H   6.924049  -2.566572  -6.208144
     C   9.425132  -4.480007  -4.820464
     H   8.099764  -4.787895  -6.494004
     H   7.613664  -1.010362  -4.372973
     C   9.165953  -2.331961  -3.603834
     C   9.834946  -3.602938  -3.770049
     H   9.940882  -5.445105  -4.937435
     C   9.575768  -1.454891  -2.553419
     C  10.887437  -3.953477  -2.870144
     C  10.604204  -1.828063  -1.692732
     H   9.060018  -0.489793  -2.436448
     C  11.265470  -3.082639  -1.852227
     H  11.387236  -4.924536  -3.000909
     H  10.901136  -1.147003  -0.879879
     H  12.076850  -3.368326  -1.165739
     C   7.735429   3.082641  -5.521656
     C   8.113463   3.953478  -4.503739
     C   8.396696   1.828065  -5.681151
     H   6.924049   3.368328  -6.208144
     H   7.613664   4.924538  -4.372974
     C   9.165953   3.602939  -3.603834
     C   9.425132   1.454893  -4.820463
     H   8.099764   1.147004  -6.494004
     C   9.834946   2.331962  -3.770049
     C   9.575768   4.480009  -2.553419
     H   9.940882   0.489795  -4.937434
     C  10.887437   1.981423  -2.870144
     C  10.604204   4.106837  -1.692732
     H   9.060018   5.445107  -2.436449
     C  11.265470   2.852261  -1.852227
     H  11.387236   1.010363  -3.000909
     H  10.901136   4.787897  -0.879879
     H  12.076850   2.566574  -1.165739
     C   1.758238   0.114256   1.835565
     C   1.383454   0.986914   0.819134
     C   1.097522  -1.140337   1.993640
     H   2.571346   0.398715   2.520618
     C   0.334430   0.635395  -0.082076
     H   1.885270   1.956398   0.685839
     C   0.069095  -1.512687   1.133501
     H   1.396710  -1.820920   2.805984
     C  -0.334354  -0.635434   0.081958
     C  -0.069069   1.512709  -1.133531
     H  -0.442641  -2.479395   1.249936
     C  -1.383484  -0.986925  -0.819111
     C  -1.097557   1.140346  -1.993567
     H   0.442653   2.479424  -1.249932
     C  -1.758283  -0.114241  -1.835508
     H  -1.885277  -1.956423  -0.685867
     H  -1.396724   1.820909  -2.805950
     H  -2.571329  -0.398703  -2.520633
     C   9.589816   0.115176   1.834689
     C   9.211776   0.986022   0.816785
     C   8.928518  -1.139375   1.994266
     H  10.401194   0.400885   2.521169
     C   8.159304   0.635481  -0.083144
     H   9.711581   1.957082   0.686036
     C   7.900153  -1.512563   1.133490
     H   9.225474  -1.820453   2.807089
     C   7.490362  -0.635521   0.083033
     C   7.749460   1.512588  -1.133509
     H   7.384382  -2.477638   1.250538
     C   6.437783  -0.986030  -0.816754
     C   6.721049   1.139380  -1.994196
     H   8.265226   2.477665  -1.250534
     C   6.059734  -0.115160  -1.834627
     H   5.938010  -1.957110  -0.686064
     H   6.424112   1.820443  -2.807055
     H   5.248421  -0.400873  -2.521182
     C   2.708486  -0.115174   5.539165
     C   3.086526  -0.986020   6.557070
     C   3.369785   1.139376   5.379589
     H   1.897108  -0.400883   4.852685
     C   4.138999  -0.635480   7.456998
     H   2.586722  -1.957081   6.687818
     C   4.398150   1.512564   6.240364
     H   3.072828   1.820455   4.566765
     C   4.807939   0.635522   7.290822
     C   4.548842  -1.512586   8.507364
     H   4.913920   2.477640   6.123316
     C   5.860520   0.986031   8.190608
     C   5.577254  -1.139379   9.368050
     H   4.033077  -2.477663   8.624388
     C   6.238568   0.115162   9.208482
     H   6.360292   1.957111   8.059918
     H   5.874190  -1.820441  10.180909
     H   7.049881   0.400875   9.895037
     C   1.032832   3.082641   9.226081
     C   1.410866   3.953478  10.243998
     C   1.694099   1.828065   9.066586
     H   0.221452   3.368328   8.539593
     H   0.911066   4.924538  10.374764
     C   2.463356   3.602939  11.143903
     C   2.722535   1.454893   9.927274
     H   1.397166   1.147004   8.253733
     C   3.132349   2.331962  10.977688
     C   2.873170   4.480009  12.194318
     H   3.238285   0.489795   9.810303
     C   4.184839   1.981423  11.877593
     C   3.901606   4.106837  13.055005
     H   2.357420   5.445107  12.311288
     C   4.562873   2.852261  12.895510
     H   4.684639   1.010363  11.746828
     H   4.198539   4.787897  13.867858
     H   5.374253   2.566574  13.581998
     C   6.059734   5.819740  -1.834628
     C   6.437783   4.948871  -0.816754
     C   6.721049   7.074280  -1.994196
     H   5.248421   5.534027  -2.521183
     C   7.490363   5.299379   0.083033
     H   5.938010   3.977790  -0.686064
     C   7.749460   7.447488  -1.133509
     H   6.424112   7.755343  -2.807055
     C   8.159304   6.570381  -0.083144
     C   7.900153   4.422337   1.133490
     H   8.265226   8.412565  -1.250534
     C   9.211776   6.920922   0.816785
     C   8.928517   4.795525   1.994266
     H   7.384382   3.457261   1.250538
     C   9.589816   6.050076   1.834689
     H   9.711581   7.891982   0.686036
     H   9.225474   4.114447   2.807089
     H  10.401194   6.335785   2.521169
     C   4.384090   3.082597   1.852325
     C   4.762144   3.953474   2.870185
     C   5.045437   1.828082   1.692675
     H   3.572778   3.368332   1.165778
     C   5.814707   3.602965   3.769996
     H   4.262367   4.924555   3.000859
     C   6.073776   1.454858   2.553450
     H   4.748475   1.147001   0.879846
     C   6.483596   2.331937   3.603858
     C   6.224526   4.480043   4.820405
     H   6.589562   0.489803   2.436348
     C   7.536158   1.981427   4.503669
     C   7.252866   4.106819   5.681180
     H   5.708740   5.445098   4.937507
     C   7.914213   2.852305   5.521529
     H   8.035936   1.010346   4.372995
     H   7.549827   4.787901   6.494008
     H   8.725524   2.566570   6.208076
     C   2.708486   5.819726   5.539165
     C   3.086526   4.948880   6.557070
     C   3.369785   7.074276   5.379588
     H   1.897108   5.534017   4.852685
     C   4.138999   5.299420   7.456998
     H   2.586722   3.977819   6.687818
     C   4.398150   7.447464   6.240364
     H   3.072828   7.755355   4.566765
     C   4.807940   6.570423   7.290822
     C   4.548842   4.422314   8.507364
     H   4.913921   8.412540   6.123316
     C   5.860520   6.920931   8.190608
     C   5.577253   4.795521   9.368051
     H   4.033077   3.457236   8.624389
     C   6.238568   6.050061   9.208482
     H   6.360293   7.892011   8.059919
     H   5.874190   4.114458  10.180909
     H   7.049881   6.335774   9.895037
     C   6.073776   7.389758   2.553450
     C   6.483596   8.266837   3.603858
     C   5.045437   7.762982   1.692674
     H   6.589563   6.424703   2.436348
     C   7.536158   7.916327   4.503669
     C   5.814706   9.537864   3.769996
     H   4.748475   7.081901   0.879846
     C   4.384090   9.017496   1.852325
     C   7.914213   8.787205   5.521529
     H   8.035936   6.945246   4.372995
     C   4.762144   9.888374   2.870185
     C   6.224526  10.414943   4.820405
     H   3.572778   9.303232   1.165778
     H   8.725524   8.501470   6.208077
     C   7.252866  10.041719   5.681180
     H   4.262367  10.859455   3.000859
     H   5.708740  11.379998   4.937507
     H   7.549827  10.722800   6.494008
     H  -0.442641   9.390404   1.249936
     C   0.069095  10.357112   1.133501
     C  -0.334354  11.234366   0.081958
     C   1.097522  10.729462   1.993640
     C  -1.383484  10.882875  -0.819111
     C   0.334430  12.505194  -0.082076
     C   1.758238  11.984056   1.835565
     H   1.396710  10.048880   2.805984
     C  -1.758283  11.755558  -1.835508
     H  -1.885277   9.913377  -0.685867
     C   1.383454  12.856714   0.819134
     C  -0.069069  13.382508  -1.133531
     H   2.571346  12.268515   2.520618
     C  -1.097557  13.010146  -1.993567
     H  -2.571329  11.471096  -2.520633
     H   1.885270  13.826197   0.685839
     H   0.442653  14.349223  -1.249932
     H  -1.396724  13.690709  -2.805950


