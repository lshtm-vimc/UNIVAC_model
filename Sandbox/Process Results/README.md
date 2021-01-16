# UNIVAC_model

## Comment by Hira Tanvir:
Some additional code to help you process the UNIVAC outcomes.
The model produces outputs by date, pathogen, coverage type and country and some other functions whether waning is turned on or off.
The code merge.psa.results.R can also be used to merge deterministic ('deter') outputs but you need to give a pattern to identify all the relevant files to merge.
i.e. if the all the files you want to merge for every country end with 'waningTRUE' then you can use this as a pattern in the function 'file_to_df'
to combine all the country outputs with that pattern identified in the output file name.
