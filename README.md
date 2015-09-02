### SNP_EM_HMM
## by Noel Carrascal
Machine Learning hack for automated single nucleotide polimorphism array data.
This software have a few bugs and works under careful selection of some of its parameters.
It is an intial proof of concept that needs to be polished and enhanced. Methodology
can be applied to other pattern recognition problems after refinement. 

Known Bugs:
* Selection of best model using vitervi's algorithm is faulty. It works most of the time but fails in some instances.

Possible Enhancements:
* Make adjustable parameter selection automatically.
* Automatically select the amount of data to be processed so that normal distributions fit the data better. In other words, data could be read in 'sliding windows' that overlap so that expectation maximization has multiple chances to fit the normal distributions. This would require some sort of concesus selection of the right fit. 
