# Analysis for 3x3

For replication on your local machine : 
- Download the `modules` folder and make sure you have an `__init__.py` file.
- Download files 1 to 6 here
- Specify the paths for data storage, `size=3` and `budget` in the `params.py` file

Then :
- *Step 1:* Run file `1_Landscapes.py` to generate all the unique landscapes and the folder structure for the storage of intermediary results
- *Step 2:* Run file `2_Optimization.py` to generate the optimal successions from the 'brute force' approach. The resulting landscapes will not necessarily be part of the initial unique landscapes and may be equivalent to other landscapes.
- *Step 3:* You may have to restart so that `params.files` gets updated with the files resulting from step 2. Then, run file `3_Keys_for_matching.py` to obtain matching keys to go from $t$ to the $t+1$ landscape, matching the succession to its equivalent.
- *Step 4:* Run `4_Successions_matched.py` to obtain the 20 period successions given a biodiversity constraint and a budget
- *Step 5:* Run `5_Statistics.py` to compute the statistics over all the datasets for a given budget
- *Step 6:* Run `6_Centrality_option.py` to compute optimal successions using the eigencentrality-based algorithm. This option is longer, as it does not rely on any matching.
As a work in progress, the `2_Optimization.py` file will soon offer a choice between methods. In this case, both 'brute_force' and 'eigencentrality based' algorithms should be usable in the routine based on matching.
