# Summary

Modification of the R package [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html) to calculate and output additional internal quantities.  Specifically $\sigma$ for each data point and the "distances" used by tSNE $p_{ij}$ and $q_{ij}$ are output.  $p_{ij}$ is returned as it is calculated (i.e., uses the Barnes-Hutt approximation to calculate $p_{ij}$).  $q_{ij}$ is always calculated exactly and in full.

The purpose of this package is mainly educational, although the additional quantities returned may be useful for interpreting points in tSNE maps more generally.  See [here]() for more details.

#Notes

The calculation isn't particularly efficient in memory or time, so probably don't use this unless you need to.

The function won't return $\sigma$ directly.  Instead it returns $beta = \frac{1}{2\sigma^2}$.
