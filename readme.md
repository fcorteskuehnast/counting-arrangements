# Requirements 

python3.8, sagemath, c++1z, gnu++1z, gmp, gmpxx, boost


# Usage

- compile the c++ code from the `cpp` directory using the command 
	```sh build.sh```
- then the main script can be ran with
	```python scripts/ipe2bound.py [config]```

where `[config]` is a configuration file such as the ones provided in the `configuration` directory.

- `constructions/k4/k4.config` encodes the patches for the $k=4$ configuration 
- `constructions/k6/k6.config` encodes the patches for the $k=6$ configuration 
- `constructions/k12rect/k12rect.config` encodes the patches for the $k=12$ configuration 

Since the computing times are several CPU days, we also decided to provide toy configurations
which can be verified with few CPU time:

- `constructions/k6/k6_small.config` encodes smaller patches for the $k=6$ configurations.

Each configuration file encodes the parameter $k$ and the areas of the regions $R_i$. 
Moreover, for each region $R_i$ it specifies a patch type $P_i$, which is encoded in an ipe-file.  
Our script reads the bipermutation from each patch type, computes the number of reroutings, 
and ultimately computes the leading constant of $B_n$ obtained by the provided configuration.


# LGV Computations

Compute the reroutings for a square patch of size $l \times l$ using the LGV lemma, 
we provide an auxiliary script `scripts/lgv.sage`.
For the config files, we provide the precomputed value for $l=500$
which can be recomputed within about 40 CPU minutes using the following command:
```sage scripts/lgv.sage 500```


# An example

```
$ python scripts/ipe2bound.py constructions/k6/k6.config 
...
============== SUMMARY ==============
        log(F(P_R4)) =                         122.94470366168639 (computing time: 18787.90s)
        F(P_R5)      =              32207077855497546508132740267 (computing time: 4219.49s)
        log(F(P_R6)) =                         102.11445434109544 (computing time: 868737.00s)
        log(F(P_R3)) =                                   349033.0 (computing time: 0.00s)
numbers of patches:
        µ_R4(P_R4, n)      = 1/4608 n² - O(n)      = 1/128 m² - O(m)   
        µ_R5(P_R5, n)      = 1/1728 n² - O(n)      = 1/48 m² - O(m)    
        µ_R6(P_R6, n)      = 1/1008 n² - O(n)      = 1/28 m² - O(m)    
        µ_R3(P_R3, n)      = 1/12000000 n² - O(n)  = 3/1000000 m² - O(m)
contributions by region:
        c_R4    ≈ 0.02668
        c_R5    ≈ 0.05480
        c_R6    ≈ 0.10130
        c_R3    ≈ 0.02909
complete bound:
        log(F_6(n))  > 0.21187484409728075
        log(B(n))    > 0.2542498129167369

============== TABLE ==============
  region    log2(F)    lim µ/n²  contribution    CPU time
----------------------------------------------------------
      R4     122.94      1/4608       0.02668   18787.90s
      R5      94.70      1/1728       0.05480    4219.49s
      R6     102.11      1/1008       0.10130  868737.00s
     R3*  349033.00  1/12000000       0.02909       0.00s
----------------------------------------------------------
       Σ          -           -       0.21187  891744.39s
```
