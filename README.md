# BinEX
Binning tool for improve binning percentage in existing binning methods

Run the script
------------

* run the python script with taxonomic label data

      $python mahalanobis_binning.py 10s binned_contigs_features.csv unbinned_contigs_features.csv taxon.csv final_ouput.csv --p 1500 9



      argument 1 : dataset name (ex : '10s')

      argument 2 : initial binning method feature file (ex : binned_contigs_features.csv)

      argument 3 : unbinned contigs feature file (ex : unbinned_contigs_features.csv)

      argument 4 : taxon file [this contain the taxons of all contigs] (ex : taxon.csv )[OPTIONAL]

      argument 5 : final output file (ex :final_ouput.csv)
      
      argument 6 : indicate parameter selcting (ex: "--p") [OPTIONAL]
      
      argument 7 : filtering length (ex: 1500) [OPTIONAL]
      
      argument 8 : Critical value (ex: 9) [OPTIONAL]


* sample dataset in shown in the test folder.

* feature files for unbinned contigs and binned contigs should follow following structure.
      
      Binned feature file headers => [contig_id, contig_length, bin , feature_1, feature_2,....., feature_N ]
      
      Unbinned feature file headers => [contig_id, contig_length, feature_1, feature_2,....., feature_N ]
      
 * taxon file contian following structure (not a requirement to have headers in csv file)
      
       Taxon file headers => [contig_id, Actual taxon]

* argument 4 (taxon file) is needed to calculate accuracy and other properties of results. Even without that it is possible to run the script and get the output.

      $python mahalanobis_binning.py 10s binned_contigs_features.csv unbinned_contigs_features.csv final_ouput.csv
      
* You can pick the filtering length and another threshold value in command line option. Depending on those values the results may changes.The order of argument follow the following.

       --p [filtering length] [critical value]
      
      eg :      
       taxon file available case :
            python mahalanobis_binning.py 10s sample_data/10s/10s_binned_contigs_features.csv  sample_data/10s/10s_unbinned_contigs_features.csv sample_data/10s/10s_taxon.csv sample_data/10s/10s_final_ouput.csv --p 1500 9
            
      taxon file not available case : 
       python mahalanobis_binning.py 10s sample_data/10s/10s_binned_contigs_features.csv  sample_data/10s/10s_unbinned_contigs_features.csv sample_data/10s/10s_final_ouput.csv --p 1500 9

* you can provide the relative paths or absolute path of files.

* final output will display following format csv file
             
      output file => [id, bin]


* PythonTestBed.java file is there to run the script using java.(usefull for web applications). If it doesn't need (for normal run), make
            
      testBedEnable = 0
      
* for degugging perposes it is possible to enbale progressbar
      
      progressBarEnable = 1
