# sQTL cluster viz

Allows for visualization of LeafCutter clusters, stratefied by genotype of a specified SNP.

sqtl_viz requires command line args:
1. Path of sqtl_viz_input_file (more on this below)
2. Path of output directory (ie where to save results).

The sqtl_viz_input_file is column-structured as follows:
1: variant id
2: cluster id
3: junction start pos
3: junction end pos
4: gene-id
5: Tissue-id
6: Path of genotype file

The sqtl_viz_input_file currently does not accept a header line.

An example of how to run sqtl_cluster_viz can be found in `driver_key.sh`.

sQTL_cluster_viz results in a pdf figure (located in the output directory) for each line in the sqtl_viz_input_file.


## Requirements
This script relies on some files (GTEx data) that is located on MARCC. And therefor cannot be run anywhere expcept on MARCC.


## Authors

* **Ben Strober** -- [BennyStrobes](https://github.com/BennyStrobes) -- bstrober3@gmail.com
