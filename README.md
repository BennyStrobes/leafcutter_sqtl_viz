# sQTL cluster viz

Allows for visualization of LeafCutter clusters, stratefied by genotype of a specified SNP.

sqtl_viz requires command line args:
1. Path of sqtl_viz_input_file (more on this below)
2. Path to output directory (ie where to save results).

The sqtl_viz_input_file is column-structured as follows:
1: variant id
2: cluster id
3: junction start pos
3: junction end pos
4: gene-id
5: Tissue-id
6: Path to genotype file

The sqtl_viz_input_file currently does not accept a header line.

An example of how to run sqtl_cluster_viz can be found in `driver_key.sh`




## Authors

* **Ben Strober** -- [BennyStrobes](https://github.com/BennyStrobes) -- bstrober3@gmail.com
