import numpy as mx
import numpy as np
import pandas as pd
from datetime import datetime
import sys
import os

def madhyper_process(prefix, outdir):
    print("start load:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    bigmas = mx.array(np.loadtxt(os.path.join(outdir, prefix+'_bigmas.tsv'), delimiter='\t', dtype=np.float32))
    bigmbs = mx.array(np.loadtxt(os.path.join(outdir, prefix+'_bigmbs.tsv'), delimiter='\t', dtype=np.float32))
    mdh = mx.array(np.loadtxt(os.path.join(outdir, prefix+'_mdh.tsv'), delimiter='\t', dtype=np.int32))
    rowinds_bigmas=mx.arange(bigmas.shape[0])
    rowinds_bigmbs=mx.arange(bigmbs.shape[0])
    results = []
    chunk_size =500  # Define your chunk size
    n_wells = bigmas.shape[1]  # Assuming n_wells is the number of columns in bigmas
    
    # Determine the total number of chunks, with a minimum of 1
    total_chunks = max(bigmas.shape[0] // chunk_size, 1)
    print('total number of chunks',total_chunks)
    b_total = mx.sum(bigmbs > 0, axis=1,keepdims=True)
    bigmbs=(bigmbs > 0).T.astype(mx.float32)
    print("start time for MAD-HYPE:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    for ch in range(0, bigmas.shape[0], chunk_size):
        print(f'Working on chunk {ch}')
        percent_complete = int((ch // chunk_size + 1) / total_chunks * 100)
        # Print progress only on 5%, 10%, etc.
        if percent_complete % 10 == 0 and percent_complete > 0:
            print(f'Progress: {ch} ({percent_complete}%)')
        chunk_end = min(ch + chunk_size, bigmas.shape[0])
        row_range = slice(ch, chunk_end)
        a_total = mx.sum(bigmas[row_range] > 0, axis=1,keepdims=True)
        overlaps = mx.matmul((bigmas[row_range] > 0).astype(mx.float32), bigmbs) #optimized
        mask_condition=-(overlaps.T - b_total).T < mdh[(overlaps).astype(mx.int16), -(overlaps - a_total).astype(mx.int16)]# mad hype only
        pairs = mx.argwhere(mask_condition)
        #pairs = mx.array(np.argwhere(np.array(mask_condition, copy=False)))
    
        result = { #this works in cupy
              'alpha_nuc': 1+(rowinds_bigmas[row_range][pairs[:, 0]]),
              'beta_nuc': 1+(rowinds_bigmbs[pairs[:, 1]]),
              'wij': (overlaps[pairs[:, 0], pairs[:, 1]]),
              'wa': (a_total[:,0][pairs[:, 0]]),
              'wb': (b_total[:,0][pairs[:, 1]])
        }

        #result = {
         #   'alpha_nuc': 1+np.array(rowinds_bigmas[row_range][pairs[:, 0]]),
          #  'beta_nuc': 1+np.array(rowinds_bigmbs[pairs[:, 1]]),
    #       'r': np.array(pairwise_cors_method2[pairs[:, 0], pairs[:, 1]]),
    #     'pval': pairwise_cors_method2_ps[pairs[:, 0], pairs[:, 1]],
    #     'pval3': pairwise_corsn[pairs[:, 0], pairs[:, 1]],
         #   'wij': np.array(overlaps[pairs[:, 0], pairs[:, 1]]),
        #    'wa': np.array(a_total[:,0][pairs[:, 0]]),
        #    'wb': np.array(b_total[:,0][pairs[:, 1]])
        #}

        results.append(result)
#result is a list of dictionaries, each dictionary contains the results for a chunk of rows. You can convert it to a pandas DataFrame like this: 
    print("end time for MAD-HYPE:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

#make pandas dataframe from each element of results and concatenate them
    results_df = pd.concat([pd.DataFrame(result) for result in results])
    results_df.to_csv(os.path.join(outdir, prefix+'_madhyperesults.csv'), index=False)
    print(f"Number of pairs: {results_df.shape[0]}")



def correlation_process(prefix, outdir, min_wells=2):
    print("start load for T-Shell:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    bigmas = mx.array(np.loadtxt(os.path.join(outdir, prefix+'_bigmas.tsv'), delimiter='\t', dtype=np.float32))
    bigmbs = mx.array(np.loadtxt(os.path.join(outdir, prefix+'_bigmbs.tsv'), delimiter='\t', dtype=np.float32))
    rowinds_bigmas = mx.arange(bigmas.shape[0])
    rowinds_bigmbs = mx.arange(bigmbs.shape[0])
    print('file read done')
    
    # Filter by min_wells
    non_zero_counts_bigmas = mx.sum(bigmas > 0, axis=1)
    non_zero_counts_bigmbs = mx.sum(bigmbs > 0, axis=1)
    valid_rows_bigmas = mx.nonzero(non_zero_counts_bigmas > min_wells)[0]
    valid_rows_bigmbs = mx.nonzero(non_zero_counts_bigmbs > min_wells)[0]
    bigmas = bigmas[valid_rows_bigmas]
    bigmbs = bigmbs[valid_rows_bigmbs]
    rowinds_bigmas = rowinds_bigmas[valid_rows_bigmas]
    rowinds_bigmbs = rowinds_bigmbs[valid_rows_bigmbs]
    
    results = []
    chunk_size = 500
    # Determine the total number of chunks, with a minimum of 1
    total_chunks = max(bigmas.shape[0] // chunk_size, 1)
    print('total number of chunks', total_chunks)
    
    # Pre-process bigmbs
    bigmb_w1_scaled = bigmbs - mx.mean(bigmbs, axis=1, keepdims=True)
    bigmb_w1_scaled = (bigmb_w1_scaled / mx.linalg.norm(bigmb_w1_scaled, ord=2, axis=1, keepdims=True)).T
    b_total = mx.sum(bigmbs > 0, axis=1, keepdims=True)
    bigmbs = (bigmbs > 0).T.astype(mx.float32)
    
    print("start processing time for T-Shell:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    for ch in range(0, bigmas.shape[0], chunk_size):
        percent_complete = int((ch // chunk_size + 1) / total_chunks * 100)
        if percent_complete % 10 == 0 and percent_complete > 0:
            print(f'Progress: {ch} ({percent_complete}%)')
            
        chunk_end = min(ch + chunk_size, bigmas.shape[0])
        row_range = slice(ch, chunk_end)
        
        # Calculate correlations
        bigma_w1_scaled = bigmas[row_range] - mx.mean(bigmas[row_range], axis=1, keepdims=True)
        bigma_w1_scaled = bigma_w1_scaled / mx.linalg.norm(bigma_w1_scaled, ord=2, axis=1, keepdims=True)
        pairwise_cors_method2 = mx.matmul(bigma_w1_scaled, bigmb_w1_scaled)
        
        # Calculate overlaps
        a_total = mx.sum(bigmas[row_range] > 0, axis=1, keepdims=True)
        overlaps = mx.matmul((bigmas[row_range] > 0).astype(mx.float32), bigmbs)
        
        # Use argsort to get top 3 correlations for each row
        sorted_indices = mx.argsort(-pairwise_cors_method2, axis=1)  # Sort in descending order
        mask = mx.zeros_like(pairwise_cors_method2, dtype=bool)
        rows, cols = mx.ogrid[:pairwise_cors_method2.shape[0], :3]
        mask[rows, sorted_indices[rows, cols]] = True
        
        # Get pairs of indices where mask is True
        pairs = mx.argwhere(mask)
        
        result = {
            'alpha_nuc': 1 + rowinds_bigmas[row_range][pairs[:, 0]],
            'beta_nuc': 1 + rowinds_bigmbs[pairs[:, 1]],
            'r': pairwise_cors_method2[pairs[:, 0], pairs[:, 1]],
            'wij': overlaps[pairs[:, 0], pairs[:, 1]],
            'wa': a_total[:,0][pairs[:, 0]],
            'wb': b_total[:,0][pairs[:, 1]]
        }
        
        results.append(result)
    
    print("end time for T-Shell:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    results_df = pd.concat([pd.DataFrame(result) for result in results])
    results_df.to_csv(os.path.join(outdir, prefix+'_corresults.csv'), index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <prefix> <outdir>")
        sys.exit(1)
        
    prefix = sys.argv[1]
    outdir = sys.argv[2]
    madhyper_process(prefix, outdir)
    correlation_process(prefix, outdir)
