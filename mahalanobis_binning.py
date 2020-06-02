import sys
import numpy as np
import pandas as pd
import scipy.linalg as sp
import matplotlib.pyplot as plt
import time as time
from scipy.stats import chi2
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib.gridspec import GridSpec
import matplotlib.cm as cm

def calculate_binwise_dist (label, df_2Tbinned) :
    arr = df_2Tbinned[df_2Tbinned['bin'] == label]
    return arr

def calculate_mahalanobis_dist (x=None, data=None, cov=None) :
    """Compute the Mahalanobis Distance between each row of x and the data  
    x    : vector or matrix of data with, say, p columns.
    data : ndarray of the distribution from which Mahalanobis distance of each observation of x is to be computed.
    cov  : covariance matrix (p x p) of the distribution. If None, will be computed from data.
    D^2 = (x-m)^T * C^(-1)* (x-m)
    """
    x_minus_mu = x - np.mean(data)
    if not cov:
        cov = np.cov(data.values.T)
    inv_covmat = sp.inv(cov)
    left_term = np.dot(x_minus_mu, inv_covmat)
    mahal = np.dot(left_term, x_minus_mu.T)
    return mahal.diagonal()

def getLableForBinsDF (binned_dataframe, taxons) :
    Allbins = binned_dataframe.drop_duplicates(subset='bin', keep="last").reset_index(drop=True)
    rows_n = Allbins['bin'].to_list()
    
    Alltaxons = taxons.drop_duplicates(subset='Actual taxon', keep="last").reset_index(drop=True)
    cols_n = Alltaxons['Actual taxon'].to_list()
    arr = np.zeros((len(rows_n),len(cols_n)))
    
    for k in range(len(rows_n)):
        for j in range(len(cols_n)):
            n = len(binned_dataframe[(binned_dataframe['bin'] == rows_n[k]) & (binned_dataframe['Actual taxon']==cols_n[j])])
            arr[k][j] =n
    
    tax = []
    for k in range(len(rows_n)):
        max = 0
        for j in range(len(cols_n)):            
            if arr[k][j] >= max :
                max = arr[k][j]
                max_col_index = j
        tax.append(cols_n[max_col_index])
      
    df_bin = pd.DataFrame(data=rows_n)
    df_bin['taxon'] = tax
    df_bin.columns = ['bin', 'taxon']
    
    return df_bin;

def run_model(argv) :

    #  EXTRACTION OF FILES & PARAMETERS
    critical_value = 9
    length_consider = 1500
    taxon_file = None
    
    if (len(argv) == 5) :
        name = argv[1]
        df_2Tbinned = pd.read_csv(argv[2])
        df_2Tunbinned = pd.read_csv(argv[3])
        output = argv[4]
        
    elif (len (sys.argv) == 6) :
        name = argv[1]
        df_2Tbinned = pd.read_csv(argv[2])
        df_2Tunbinned = pd.read_csv(argv[3])
        taxon_file = pd.read_csv(argv[4], names=['id','Actual taxon'], header=None)
        output = argv[5]
    else:        
        print ("Incorrect arguments / Missing arguments")
        exit()
    
    # COMPUTATION

    df= df_2Tunbinned # for debugging   
    bins = df_2Tbinned['bin'].unique() # get the bins
    bins_array = []
    ignore_bins = []

    for bin in bins : # preposessing part handle
        if ((df_2Tbinned['bin'] == bin).sum()> 20):
            bins_array.append(bin)
        else:
            ignore_bins.append(bin)

    i = 0
    newdf = pd.DataFrame(columns=['id', 'ofdeg', 'gc', 'length','bin']) # this df contain newly binned contig
    mahala_dist = []

    while i < len(df.index):
        row = df.loc[i,]
        lowest_variance = sys.float_info.max
        assigned_bin = None
        goodForBin = 0
        distance = None
        minD = sys.float_info.max; # this variable for store each unbinned contigs mahalanobis dist
        
        for bin in bins_array :
            label_bin = calculate_binwise_dist(bin, df_2Tbinned)
            df1 = label_bin[['id','length', 'ofdeg', 'gc','principal component 1', 'principal component 2','principal component 3']].append(row, ignore_index=True) # add each unbinned contigs to bin and calculate distance
            df2 = df1[['ofdeg', 'gc','principal component 1', 'principal component 2','principal component 3']]
            df1['mahala'] = calculate_mahalanobis_dist(x=df2, data=df2[['ofdeg', 'gc','principal component 1', 'principal component 2','principal component 3']]) #dataframe with distance column
            dist = df1.loc[df1.index[-1], "mahala"]
            d = dist
            if d < minD : # get the minimum mahalanobis dist
                minD = d
            
            if dist < critical_value : # this will check whether contig has smaller distance or not
                goodForBin = 1
                x = dist; # contig mahalanobis distances
                m = df1['mahala'].mean();  # mean
                n = len(df1['mahala'].index); #sample number 
                variance = (x-m)*(x-m) / n 
                
                if variance < lowest_variance :
                    lowest_variance = variance
                    assigned_bin = bin

    if newdf.empty:
        print('Does not bin any contig')

    else :
        final_out = pd.concat([newdf[['id','bin']], df_2Tbinned[['id','bin']]], ignore_index=True)
        final_out.to_csv(output, index=False,header=True)
        
    if taxon_file is None :
        print ("finish binning")
        exit()
    else :        
    
        binned_data_with_taxon = df_2Tbinned.merge(taxon_file, on="id", how = 'left')
        bins_from_intial = getLableForBinsDF(binned_data_with_taxon,taxon_file) #rows bins, columns taxons

        df_2Tbinned = df_2Tbinned.merge(bins_from_intial, on="bin", how = 'inner')
        df_2Tbinned = df_2Tbinned.merge(taxon_file, on="id", how = 'inner')
        df_2Tbinned['Is prediction correct']= (df_2Tbinned['taxon']==df_2Tbinned['Actual taxon'])

        true_prediction_2t = df_2Tbinned['Is prediction correct'].values.sum() # true count
        false_prediction_2t = (~df_2Tbinned['Is prediction correct']).values.sum() # false count
        binned_count_2t = true_prediction_2t + false_prediction_2t

        accuracy_2t = 100 * true_prediction_2t / binned_count_2t
        print("Binning accuracy in 2t method:% 5.2f" % accuracy_2t, "%")

        per_of_contig_binned_2t = 100 * len(df_2Tbinned.index)/(len(df_2Tunbinned.index) +len(df_2Tbinned.index))
        print("percentage contig bin in 2t :% 5.2f" % per_of_contig_binned_2t,"%")
        
        if newdf.empty:
            print('Does not bin any contigs in our method')

        else :        
            
            # check if file is given    
            results = newdf.merge(bins_from_intial, on="bin", how = 'left').reset_index(drop=True)
            results = results.merge(taxon_file, on="id", how = 'left').reset_index(drop=True)

            results = results[['id','length','bin','mahala','taxon','Actual taxon']]
            results['Is prediction correct']= (results['taxon']==results['Actual taxon'])

            results = results[['id','length','bin','mahala','taxon','Actual taxon']]
            results['Is prediction correct']= (results['taxon']==results['Actual taxon'])

            true_prediction = results['Is prediction correct'].values.sum() # true count
            false_prediction = (~results['Is prediction correct']).values.sum() # false count
            binned_count_our_model = false_prediction + true_prediction

            accuracy = 100 * (true_prediction + true_prediction_2t) / (binned_count_our_model + binned_count_2t)
            print("Binning accuracy in proposed model:% 5.2f" % accuracy, "%")

            per_of_contig_binned = 100 * (len(newdf.index)+len(df_2Tbinned.index))/(len(df_2Tunbinned.index) +len(df_2Tbinned.index))
            print("percentage contig bin in 2T and our model :% 5.2f" % per_of_contig_binned,"%")
        
              
# single running run the script as following
# python mahalanobis_binning.py 10s sample_data/10s/10s_binned_contigs_features.csv  sample_data/10s/10s_unbinned_contigs_features.csv sample_data/10s/10s_taxon.csv sample_data/10s/10s_final_ouput.csv

if __name__ == '__main__':
    start = time.time()
    run_model(sys.argv)        
    end_time = time.time()
    t = start-end
    print("Average time taken for execution "+ str(int(t/60))+" min " +str(int(t%60))+" sec")


