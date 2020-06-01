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

def run_model(name,arg1,arg2,arg3,arg4,arg5,arg6,arg7,out_file) :

    # Assuming it follows chi square distribution with 0.01 significance level and 4 degrees of freedom
    # Critical values for two degrees of freedom

    critical_value = chi2.ppf((1-0.01), df=2)
    print(critical_value);
    
    #  EXTRACTION OF FILES & PARAMETERS
    critical_value = 9
    length_consider = 1500
    name = '10s'
    name1 = '10s_new'

    df_2Tunbinned = pd.read_csv('sample_data/'+name1+'/'+name+'_unbinned_contigs_features.csv') # unbinned file
    df_2Tbinned = pd.read_csv('sample_data/'+name1+'/'+name+'_binned_contigs_features.csv') # binned file
    taxon_file = pd.read_csv('sample_data/'+name1+'/'+name+'_taxon.csv', names=['id','Actual taxon'], header=None) # unbinned file
   
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
            label_bin = calculate_binwise_dist(bin)
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
        print('Does not bin any contig')

    else :
        final_out = pd.concat([newdf[['id','bin']], df_2Tbinned[['id','bin']]], ignore_index=True)
        final_out.to_csv('sample_data/'+name1+'/'+name+'_final_ouput.csv', index=False,header=True)
        
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
        
        
        
#single running

arg1 ='sample_data/cami2/cami2_unbinned_contigs.OFDEG'
arg2 ='sample_data/cami2/cami2_unbinned_contigs.n4'
arg3 ='sample_data/cami2/binned_points.csv'
arg4 ='sample_data/cami2/cami2_view2.n4'
arg5 ='sample_data/cami2/sim.contig.ans'
arg6 ='sample_data/cami2/taxon_bins.csv'
arg7 ='sample_data/cami2/cami2_output.L2_BINS'

arg8 ='sample_data/cami2/cami2_model_results.csv'

# run the script as following
#python mahalanobis_binning.py 10s sample_data/10s/10s_unbinned_contigs.OFDEG sample_data/10s/10s_unbinned_contigs.n4 sample_data/10s/binned_points.csv sample_data/10s/10s_view2.n4 sample_data/10s/sim.contig.ans sample_data/10s/taxon_bins.csv sample_data/10s/10s_output.L2_BINS sample_data/10s/10s_model_results.csv

start = time.time()
# run_model(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9]) 
run_model('cami2',arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8) 
end = time.time()
t = start-end
print("Average time taken for execution "+ str(int(t/60))+" min " +str(int(t%60))+" sec")




## for bin comparison ,'sharon','25s','cami2','10s'
#dataSets = ['cami','sharon','25s','cami2','simBG','10s']
#x = np.arange(len(dataSets))
#ys = [i+x+(i*x)**2 for i in range(len(dataSets))]
#colors = cm.rainbow(np.linspace(0, 1, len(dataSets)))
#c = 0
#for name in dataSets : 
#    #ret_df = run_model(ds) # pass the dataset to model and get result dataframe for analyse
#    
#    binnedFile = pd.read_csv("sample_data/"+name+"/binned_points.csv") # binned file
#    df_2Tbinned = binnedFile[['id', 'ofdeg', 'gc', 'length','bin']]
#    file4 = pd.read_csv('sample_data/'+name+'/'+name+'_view2.n4', delimiter = ',') # file contain details about all the contigs
#    df_tnf_binned = file4.drop(['length', 'ns'], axis=1)
#
#    # tnf pca for binned data
#    ids = df_tnf_binned.loc[:,['id']] # Separating out the ids
#    tnt_frquencies = df_tnf_binned.drop(['id'], axis=1).values # Separating out the tnf freq
#    tranformed_tnf = StandardScaler().fit_transform(tnt_frquencies) # Standardizing the features
#    pca = PCA(n_components=3)
#    principalComponents = pca.fit_transform(tranformed_tnf)
#    principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2','principal component 3'])
#    principalDf['id'] = ids
#    df_2Tbinned = df_2Tbinned.merge(principalDf, on="id", how = 'inner') # add principle component cols
#    ret_df = df_2Tbinned   
#    
#    df = ret_df[['ofdeg', 'gc','principal component 1', 'principal component 2','principal component 3']]
#    ret_df['mahala'] = calculate_mahalanobis_dist(x=df, data=df[['ofdeg', 'gc','principal component 1', 'principal component 2','principal component 3']])
#    d_df = ret_df.groupby(['bin']).mean()
#    d_df = d_df[['mahala']]
#    plt.plot(d_df.index,  d_df['mahala'], color=np.random.rand(3,) , marker='o')
#    c=+1
#
#plt.xlabel('Bins')
#plt.ylabel('Distanes')
#plt.title('Distance analyse by bin')
#plt.legend(dataSets, loc='upper left')
#plt.savefig("bin_analyse.png")
#plt.show()

# # for time comparison
# n=1 #number of iteractions
# time_comp = open('time.txt', 'w')
# dataSets = ['10s','25s']
# for ds in dataSets : 
    # sum =0
    # for i in range(n):
        # start = time.time()
        # ret_df = run_model(ds) # pass the dataset to model and get result dataframe for analyse
        # end = time.time()
        # sum = sum +(start-end)

    # t_avg = sum / n 
    # time_comp.write("Average time taken for execution for " +ds+ " dataset "+ str(int(t_avg/60))+" min " +str(int(t_avg%60))+" sec/n")
    # print("Average time taken for execution "+ str(int(t_avg/60))+" min " +str(int(t_avg%60))+" sec")
    
# # CHECK HOW CLOSE THE REFERENCE BINS ARE
# # print(df_analyse)
# df2 = df_2Tbinned[['ofdeg', 'gc','principal component 1', 'principal component 2','principal component 3']]
# df_2Tbinned['mahala'] = calculate_mahalanobis_dist(x=df2, data=df2[['ofdeg', 'gc','principal component 1', 'principal component 2','principal component 3']])
# d_df = df_2Tbinned[['id','bin','mahala']]
# #print(d_df)
# def pad_or_truncate(some_list, target_len):
    # return some_list[:target_len] + [0]*(target_len - len(some_list))

# cols = d_df['bin'].unique()
# cols.sort()
# df_bins_dist = pd.DataFrame(columns=cols)
# max=0
# for col in cols:
    # x = d_df.loc[d_df['bin'] == col]['mahala'].tolist()    
    # if (len(x) > max):
        # max = len(x)

# for col in cols:
    # x = pad_or_truncate(d_df.loc[d_df['bin'] == col]['mahala'].tolist(),max)
    # df_bins_dist[col] = x
# df_bins_dist.replace(0, np.nan, inplace=True)
# ax = df_bins_dist.plot.box(grid='True',figsize=(16,8),fontsize=13, color='r')  
# ax.set_ylabel('Mahalanobis distance value')
# ax.set_xlabel('Bins')
# ax.set_title('Contig Assignment for '+name+ ' Dataset')
# #plt.savefig('3.jpg')