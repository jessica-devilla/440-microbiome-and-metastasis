data(throat.otu.tab)
data(throat.tree)
data(throat.meta)

comm <- t(throat.otu.tab)
meta.dat <- throat.meta

set.seed(123)
# For count data
zico.obj <- ZicoSeq(meta.dat = meta.dat, feature.dat = comm, 
                    grp.name = 'SmokingStatus', adj.name = 'Sex', feature.dat.type = "count",
                    # Filter to remove rare taxa
                    prev.filter = 0.2, mean.abund.filter = 0,  max.abund.filter = 0.002, min.prop = 0, 
                    # Winsorization to replace outliers
                    is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                    # Posterior sampling to impute zeros
                    is.post.sample = TRUE, post.sample.no = 25, 
                    # Multiple link functions to capture diverse taxon-covariate relation
                    link.func = list(function (x) x^0.25, function (x) x^0.5, function (x) x^0.75), 
                    stats.combine.func = max,
                    # Permutation-based multiple testing correction
                    perm.no = 99,  strata = NULL, 
                    # Reference-based multiple stage normalization
                    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                    # Family-wise error rate control
                    is.fwer = FALSE,
                    verbose = TRUE, return.feature.dat = FALSE)

