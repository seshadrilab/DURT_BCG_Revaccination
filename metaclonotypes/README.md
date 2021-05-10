# James_et_al_2021
Scripts for Durable expansion of TCR-δ meta-clonotypes after BCG revaccination in humans

### Files
```
|   (R Scripts Run Using R)
├── 2021_01_13_delta_BCG_beta_binomial_results.tsv -- regression results
├── 2021_01_13_delta_BCG_metaclones.tsv            -- definition of meta-clonotypes generated from expanded clones
├── 2021_01_13_delta_BCG_metaclones.tsv.tab.tsv    -- tabulation of TCRs conformant to meta-clonotypes in bulk files
├── README.md
├── beta_binomial_regressions.R                    -- script to run regressions from (2021_01_13_delta_BCG_metaclones.tsv.tab.tsv)
└── python (Run in Python 3.8.5 [Clang 10.0.0 ])
    ├── metaclonotypes.py    -- discover meta-clonotypes using tcrdist3
    ├── preprocess._bulk.py  -- preprocess Adaptive files for tcrdist3
    └── tabulate.py          -- tabulate TCRs conformant to meta-clonotypes in bulk files
```

##### Python Version

```
Python 3.8.5 (default, Sep  4 2020, 02:22:02)
[Clang 10.0.0 ] :: Anaconda, Inc. on darwin
Type "help", "copyright", "credits" or "license" for more information.
```

##### R Version
```
platform       x86_64-apple-darwin15.6.0   
arch           x86_64                      
os             darwin15.6.0                
system         x86_64, darwin15.6.0        
status                                     
major          3                           
minor          6.0                         
year           2019                        
month          04                          
day            26                          
svn rev        76424                       
language       R                           
version.string R version 3.6.0 (2019-04-26)
nickname       Planting of a Tree    
``