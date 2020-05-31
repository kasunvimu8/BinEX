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

def run_model(name,arg1,arg2,arg3,arg4,arg5,arg6,arg7,out_file) :

    # Assuming it follows chi square distribution with 0.01 significance level and 4 degrees of freedom
    # Critical values for two degrees of freedom

    critical_value = chi2.ppf((1-0.01), df=2)
    print(critical_value);

    # UNBINNED DATA FEATURE EXTRACTION
    critical_value = 9
    length_consider = 1500
    # name = 'sharon'

    contig_nameread = ['id', 'Actual taxon','option','len']
    file1 = pd.read_csv(arg1, delimiter = ',') # file contain details about all the contigs (argument 1)
    unbinned_contig_count = len(file1.index)
    dropindexNames1 = file1[ file1['length'] < length_consider ].index
    file1.drop(dropindexNames1 , inplace=True)
    file1.reset_index(inplace = True, drop = True)
    # print(file1)

    file2 = pd.read_csv(arg2, delimiter = ',') # file contain details about all the contigs (argument  2)
    dropindexNames2 = file2[ file2['length'] < length_consider ].index
    file2.drop(dropindexNames2 , inplace=True)
    file2.reset_index(inplace = True, drop = True)
    # print(file2)
    df_tnf = file2.drop(['length', 'ns'], axis=1)
    # print(df_tnf)

    df_2Tunbinned = file1[['id', 'ofdeg', 'gc', 'length']]

    # tnf pca for unbinned data
    ids = df_tnf.loc[:,['id']] # Separating out the ids
    tnt_frquencies = df_tnf.drop(['id'], axis=1).values # Separating out the tnf freq
    tranformed_tnf = StandardScaler().fit_transform(tnt_frquencies) # Standardizing the features
    pca = PCA(n_components=3)
    principalComponents = pca.fit_transform(tranformed_tnf)
    principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2','principal component 3'])
    principalDf['id'] = ids
    df_2Tunbinned = df_2Tunbinned.merge(principalDf, on="id", how = 'inner') # add principle component cols
    # df_2Tunbinned = file1[['id', 'ofdeg', 'gc', 'length']];
    # print(df_2Tunbinned.head(10))

    # BINNED DATA FEATURE EXTRACTION

    binnedFile = pd.read_csv(arg3) # binned file (argument 3)
    df_2Tbinned = binnedFile[['id', 'ofdeg', 'gc', 'length','bin']]
    file4 = pd.read_csv(arg4, delimiter = ',') # file contain details about all the contigs (argument 4)
    df_tnf_binned = file4.drop(['length', 'ns'], axis=1)

    # tnf pca for binned data
    ids = df_tnf_binned.loc[:,['id']] # Separating out the ids
    tnt_frquencies = df_tnf_binned.drop(['id'], axis=1).values # Separating out the tnf freq
    tranformed_tnf = StandardScaler().fit_transform(tnt_frquencies) # Standardizing the features
    pca = PCA(n_components=3)
    principalComponents = pca.fit_transform(tranformed_tnf)
    principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2','principal component 3'])
    principalDf['id'] = ids
    df_2Tbinned = df_2Tbinned.merge(principalDf, on="id", how = 'inner') # add principle component cols
    # print(df_2Tbinned.head(10))

    file3 =  pd.read_csv(arg5, delimiter = '\t', names=contig_nameread, header=None) # file contain names (argument 5)
    df_names = file3.iloc[: , [0, 1]]

    binned_contig_count = len(binnedFile.index)

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
    
    #print(len(ignore_bins))
    #print(ignore_bins)
    i = 0
    newdf = pd.DataFrame(columns=['id', 'ofdeg', 'gc', 'length','bin'])
    mahala_dist = []

    while i < len(df.index): 
         row = df.loc[i,]
         lowest_variance = sys.float_info.max
         assigned_bin = None
         goodForBin = 0
         distance = None
         minD = sys.float_info.max # this variable for store each unbinned contigs mahalanobis dist
        
         for bin in bins_array :
             label_bin = calculate_binwise_dist(bin,df_2Tbinned)
             df1 = label_bin[['id', 'ofdeg', 'gc', 'length','principal component 1', 'principal component 2','principal component 3']].append(row, ignore_index=True) # add each unbinned contigs to bin and calculate distance
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
                     distance = dist
         i += 1
        
         if goodForBin == 1 : # contigs bin in our model
             row['bin'] = assigned_bin
             row['mahala'] = distance
             newdf = newdf.append(row, ignore_index=True)
        
         mahala_dist.append(minD)

    if newdf.empty:
        print('Does not bin any contig')
    
    else :
        newdf[['id','bin']].to_csv(out_file, index=False)
        newdf = newdf.merge(df_names, on="id", how = 'left')

        # print(newdf)
        # print(newdf[['id', 'bin','mahala', 'Actual taxon']]); 
        pd.options.mode.chained_assignment = None  # default='warn'
        df['mahala'] = mahala_dist # assign distance for to column 'mahala' in df dataframe
        # print(df)
        # print(newdf)

        # check results for dataset

        summery = binnedFile.drop_duplicates(subset='bin', keep="last")
        summery_of_bins = summery[['bin','taxon']]
        summery_of_bins.reset_index(inplace = True, drop = True)
        # print(summery_of_bins)
        # print(type(summery_of_bins.iloc[0,0]))

        # get the bin name as most number of contig reference
        fname = 'sample_data/'+name+'/taxon_bins.csv'
        # data = pd.read_csv(fname)                            
        data = pd.read_csv(fname, index_col=0)
        bins = list(data.columns)

        tax = []
        for bin in bins:
            tax.append(data[bin].idxmax())

        df_bin = pd.DataFrame(data=bins)
        df_bin['taxon'] = tax
        df_bin.columns = ['bin', 'taxon']
        df_bin['bin'].astype('int64')
        # print (df_bin) 
        # print(type(summery_of_bins.iloc[0,0]))

        # there is problem, should solve

        # print(df_bin)
        # print(summery_of_bins)
        # print(df_bin)
        # print(newdf)
        df_bin['bin'] = df_bin['bin'].astype(int)
        results = newdf.merge(df_bin, on="bin", how = 'inner')

        results = results[['id','length','bin','mahala','taxon','Actual taxon']]
        results.columns = ['id','length','bin','mahala','taxon','Actual taxon']
        results['Is prediction correct']= (results['taxon']==results['Actual taxon'])


        results = results[['id','length','bin','mahala','taxon','Actual taxon']]
        results['Is prediction correct']= (results['taxon']==results['Actual taxon'])
        # print(results)

        true_prediction = results['Is prediction correct'].values.sum() # true count
        false_prediction = (~results['Is prediction correct']).values.sum() # false count
        binned_count_our_model = false_prediction + true_prediction
        # print(binned_count_our_model)

               
    results_2t_df = pd.read_csv('sample_data/'+name+'/'+name+'_output.L2_BINS', delimiter = '\t', names=['id', 'bin'])
    # print(results_2t_df)
    results_2t_df = results_2t_df.merge(df_bin, on="bin", how = 'inner')
    results_2t_df = results_2t_df.merge(df_names, on="id", how = 'inner')
    results_2t_df['Is prediction correct']= (results_2t_df['taxon']==results_2t_df['Actual taxon'])
    # print(results_2t_df)
    true_prediction_2t = results_2t_df['Is prediction correct'].values.sum() # true count
    false_prediction_2t = (~results_2t_df['Is prediction correct']).values.sum() # false count
    binned_count_our_model_2t = true_prediction_2t + false_prediction_2t
    
    print(results_2t_df)
    accuracy_2t = 100 * true_prediction_2t / binned_count_our_model_2t
    print("Binning accuracy in 2t method:% 5.2f" % accuracy_2t, "%")

    ac1 = 100 * len(df_2Tbinned.index)/(len(df_2Tunbinned.index) +len(df_2Tbinned.index))
    print("percentage contig bin in 2T:% 5.2f" % ac1,"%")
    
    if (newdf.empty) :
        print("Does not bin any contig.Accuracy and % bin are unchanged")
    else :
        print(results)
        ac1 = 100 * (len(newdf.index)+len(df_2Tbinned.index))/(len(df_2Tunbinned.index) +len(df_2Tbinned.index))
        print("percentage contig bin in 2T and our model :% 5.2f" % ac1,"%")
        
        accuracy = 100 * (true_prediction+true_prediction_2t) / (binned_count_our_model+binned_count_our_model_2t)
        print("Binning accuracy in our model:% 5.2f" % accuracy, "%")
    
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