"""
Program for Mouysset, Jean, 2023 - Landscape Connectivity Dilemma

STEP 1: GENERATE ALL UNIQUE LANDSCAPES IN 3X3

Requirement:
- specify path to save data in params.py
- specify path to save future results in params.py

Description:
- 1. Generate all folders for saving data and results
- 2. There are 3**9 potential landscapes. However, these can be reduced by considering :
    - Vertical and horizontal symmetries
    - 90° rotations

Outputs :
- Folder architecture for project
- All landscapes in the sense there are no symmetric or rotation equivalent in the data
"""

#region Import packages and functions
import time
import multiprocessing as mp
import pickle

#endregion

#region Import modules
exec(open('current_scripts/params.py').read())
import modules.params as params

import modules.utilities as utilities
#endregion
#region Get unique landscapes as dictionaries and save them

__spec__=None
__name__="__main__"
debut = time.time()

n_process = 14 #number of processors for multiprocessing

if __name__=="__main__":
    size=9
    list_non_z = list(range(size+1))
    list_non_z.remove(0)
    store = dict()


    #Loop for number of non-zero elements in the landscape
    for i in list_non_z:
        all_two = list(range(i+1))
        # Number of possible 2 (including 0)
        non_zz    = [i]*(i+1)
        # Number of non_zero elements
        start_time = time.time()

        # Multiprocessing : given number of non zero elements, compute all unique landscape with twos ranging from 0 to the number of non zero elements
        # which are 1s otherwise.
        with mp.Pool(processes=n_process) as pool:
            results = pool.starmap(utilities.low_lev_land3,list(zip(all_two,non_zz)))
        pool.terminate()
        print("Step :",i)
        print("took %s seconds" % (time.time() - start_time))

        # Saving output:
        a_file = open(params.path_to_save_data + '/land3_' + str(i) + ".pkl", "wb")
        pickle.dump(results, a_file)
        a_file.close()

print("Overall time : %s",(time.time()-debut))

#endregion