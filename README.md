# kb_Amplicon

This is a [KBase](https://kbase.us) module generated by the [KBase Software Development Kit (SDK)](https://github.com/kbase/kb_sdk).

This module intends to wrap the R libraries Vegan, Paramigene and so on with python.

From library R Vegan we provide the metaMDS for analyzing ecological data.

The metaMDS routine in R vegan library has the useful default behavior of following the ordination with a rotation via principal components analysis such that MDS axis 1 reflects the principal source of variation, and so on, as is characteristic of eigenvalue methods.

Procrustes rotation--The minimum sum of squares is called the root mean square error (rmse). The smaller the rmse, the more similar the two configurations are. Two final configurations are considered to have converged (arrived at essentially the same solution) when the rmse is less than 0.01, and no single residual value exceeds 0.005. Procrustes analysis thereby provides a mechanism for determining when to stop repeatedly re-running the analysis - stop when there is convergence as measured by procrustes rmse.

For ecological data, samples should be standardized by sample size to avoid ordinations that reflect primarily sample size unless the input matrix has already been standardized.  The Wisconsin double standarization is done automatically following the square transformation. 

Generate a distance (dissimiliarity) matrix from the multivariate data:
if distance_metric == 'bray', call metaMDS with the default Bray distance setting (that generates the Bray-Curtis dissimilarity matrix)
else call metaMDS with the euclidean distance matrix


You will need to have the SDK installed to use this module. [Learn more about the SDK and how to use it](https://kbase.github.io/kb_sdk_docs/).

You can also learn more about the apps implemented in this module from its [catalog page](https://narrative.kbase.us/#catalog/modules/kb_Vegan) or its [spec file]($module_name.spec).

# Setup and test

Add your KBase developer token to `test_local/test.cfg` and run the following:

```bash
$ make
$ kb-sdk test
```

After making any additional changes to this repo, run `kb-sdk test` again to verify that everything still works.

# Installation from another module

To use this code in another SDK module, call `kb-sdk install kb_Vegan` in the other module's root directory.

# Help

You may find the answers to your questions in our [FAQ](https://kbase.github.io/kb_sdk_docs/references/questions_and_answers.html) or [Troubleshooting Guide](https://kbase.github.io/kb_sdk_docs/references/troubleshooting.html).
