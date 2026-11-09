# Analysis of results

Analysis notebooks and supporting code to analyze the results coming out of the experimental nextflow pipeline described above.

## Data processing
1. Run `standardize_heatmap_data.py`. Use `--help` for usage. This takes the variable heatmap formats from all the different tools and transforms and saves them into a standardized format. The data is not scaled in any way, nor is missing data imputed.
2. Run `standardize_sample_predictions.py`. Use `--help` for usage. Depending on the tool, this takes either sample match predictions and transforms them into a standardized format or it takes the standardized heatmap data and applies a threshold to create predictions and then store that in a standardized format (not all tools output predictions). The prediction matrices label matches as 1, mismatches as 0, and unknown results as NaN.
3. Run `accuracy_calculations.py`. Use `--help` for usage. This calculates a few performance metrics such as precision, recall, and F1 score for pseudobulk experiments, and stores the results.

### Testing
Each of the scripts described above has a corresponding script of tests:
- `test_standardize_heatmap_data.py`
- `test_standardize_sample_predictions.py`
- `test_accuracy_calculations.py`

These use some fake test data (bespoke for each tool) found in `fixtures/standardize_heatmap`. 

Run the tests from this directory with `python -m unittest discover -s tests -v`.


## Analysis notebooks
- **pseudobulk_results.ipynb:** analyzes the results from the basic pseudobulk experiments. This focuses on the accuracy of the tools that we're comparing.
- **pathological_pseudobulk_experiments.ipynb:** analyzes the results from the pathological pseudobulk experiments (missing samples, double samples, uneven pseudobulk sizes, pseudobulk vs single-cell, and HYSYS heatmap + Vireo algorithm) to sanity check our results and test the limits of the tools. 
- **real_data_results.ipynb:** analyzes results from real data experiments. In any given experiment two modalities are compared and we want to check if the sample ids from those two modalities point to samples from the same patient. If not, this indicates there was a sample swap or mislabelling. 
- **computational_burden.ipynb:** compiles the computational burden results from the HGSOC (new) real data experiment using the trace files from the Nextflow runs.
- **dendrogram.ipynb:** creates dendrograms based on either the sample similarity heatmaps or the sample match predictions to assess how similar the results are from different tools.
- **matrix_dot_plot.py:** creates a matrix dot plot for all benchmarked tools from an array of values. 


## Supporting code
- **analysis_shared_functions.py:** has mostly data loading and manipulation stuff functions.
- **visualization.py:** has functions that create figures.
- **parameters.py:** contains shared parameters.

A conda environment for everything is found in **analysis.yml**.
