# Datadriven_Metagenomic_binning
Binning tool for improve binning percentage in existing binning methods

Run the script
------------

* run the python script without taxonomic label data

      $python mahalanobis_binning.py 10s binned_contigs_features.csv unbinned_contigs_features.csv taxon.csv final_ouput.csv



      argument 1 : dataset name (ex : '10s')

      argument 2 : initial binning method feature file (ex : binned_contigs_features.csv)

      argument 3 : unbinned contigs feature file (ex : unbinned_contigs_features.csv)

      argument 4 : taxon file (ex : taxon.csv)

      argument 5 : final output file (ex :final_ouput.csv)


* sample dataset in shown in the test folder.

* feature files for unbinned contigs and binned contigs should follow following structure.
      
      For binned feature file headers => [contig_id, contig_length, bin , feature_1, feature_2,....., feature_N ]
      
      For unbinned feature file headers => [contig_id, contig_length, feature_1, feature_2,....., feature_N ]

* argument 4 (taxon file) is needed to calculate accuracy and other properties of results. Even without that it is possible to run the script and get the output.

      $python mahalanobis_binning.py 10s binned_contigs_features.csv unbinned_contigs_features.csv final_ouput.csv

* you can provide the relative paths or absolute path of files.
