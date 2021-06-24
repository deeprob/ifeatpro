# ifeatpro (Physicochemical Feature Encoder for Protein Sequences)
A python package that generates 21 numerically encoded feature representation for protein sequences based on their 
physicochemical properties. 

**Note: ifeatpro is based on iFeature, a python based toolkit available at 
[link](https://github.com/Superzchen/iFeature/). Here, we have packaged 21 alignment free feature encoding functions 
available in iFeature into a pip installable module for easy usage and improved accessibility of a protein feature 
encoding tool.** 

# ifeatpro installation
```python
pip install ifeatpro
```


# ifeatpro usage
```python
from ifeatpro.features import get_feature, get_all_features
```

## Generating some random protein sequences and storing them in fasta format


```python
import random


AA = "ACDEFGHIKLMNPQRSTVWY"

sequences = ["".join([random.choice(AA) for _ in range(150)]) for _ in range(5)]

!mkdir -p ifeatpro_data

fasta_file = "ifeatpro_data/seq.fa"
with open(fasta_file, 'w') as f:
    for i, seq in enumerate(sequences):
        f.write(f">enz_{i}")
        f.write("\n")
        f.write(seq)
        f.write("\n")
```

## Getting all 21 feature encodings from protein sequences using ifeatpro 

ifeatpro contains 21 features which are capable of numerically encoding protein sequences based on their physicochemical
 properties. They are:

1. aac
2. apaac
3. cksaagp
4. cksaap
5. ctdc
6. ctdd
7. ctdt
8. ctriad
9. dde
10. dpc
11. gaac
12. gdpc
13. geary
14. gtpc
15. ksctriad
16. moran
17. nmbroto
18. paac
19. qsorder
20. socnumber
21. tpc

Using *get_all_features* function, an user can create all the 21 physicochemical encoding based feature extraction techniques provided by ifeatpro. The first argument of this function denotes the fasta file that contains protein sequences while the second argument denotes the output directory where the files will be stored as csv files. 


```python
get_all_features(fasta_file, "./ifeatpro_data/")
```


## Creating a single feature encoding using ifeatpro
An user can also create any one of the 21 feature extraction techniques available in ifeatpro using the 
*get_feature* function. The function takes the fasta file as the first argument, feature encoding type as the 
second argument and output directory where the file will be stored as the third argument. For example if an user 
wants to create aac type feature encoding using the fasta_file that we created above and would like to store it in 
ifeatpro_data directory, they can run the following command:


```python
get_feature(fasta_file, "aac", "ifeatpro_data/")
```

# feature extraction techniques description

To get a detailed description of the feature extraction techniques used in ifeatpro, please refer to the 
Supplementary Document of the paper [link to be added soon]().  

# similar modules to encode protein sequences
Other modules that can be used to generate numerical encoding of protein sequences are:
1. ngrampro [link](https://pypi.org/project/ngrampro/)
2. pssmpro [link](https://pypi.org/project/pssmpro/)

