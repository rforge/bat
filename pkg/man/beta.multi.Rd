% Generated by roxygen2 (4.0.1): do not edit by hand
\name{beta.multi}
\alias{beta.multi}
\title{Beta diversity among multiple sites.}
\usage{
beta.multi(comm, tree, abund = FALSE, func = "jaccard", raref = 0,
  runs = 1000)
}
\arguments{
\item{comm}{A sites x species matrix, with either abundance or incidence data.}

\item{tree}{An hclust or phylo object (used only for PD or FD).}

\item{abund}{A boolean (T/F) indicating whether abundance data should be used or converted to incidence before analysis.  If not specified, default is FALSE.}

\item{func}{Indicates whether the Jaccard or Soerensen family of beta diversity measures should be used. If not specified, default is jaccard.}

\item{raref}{An integer specifying the number of individuals for rarefaction (individual based).
If raref < 1 no rarefaction is made.
If raref = 1 rarefaction is made by the minimum abundance among all sites.
If raref > 1 rarefaction is made by the abundance indicated.
If not specified, default is 0.}

\item{runs}{Number of resampling runs for rarefaction. If not specified, default is 1000.}
}
\value{
A matrix of beta measures x diversity values (average and variance).
}
\description{
Beta diversity with possible rarefaction - multiple sites measure calculated as the average or variance of all pairwise values.
}
\details{
Beta diversity of multiple sites simultaneously is calculated as either the average or the variance among all pairwise comparisons (Legendre, 2014).
The beta diversity measures used here follow the partitioning framework independently developed by Podani & Schmera (2011) and Carvalho et al. (2012)
and later expanded to PD and FD by Cardoso et al. (2014), where Btotal = Brepl + Brich. Btotal = total beta diversity, reflecting both species replacement and loss/gain;
Brepl = beta diversity explained by replacement of species alone; Brich = beta diversity explained by species loss/gain (richness differences) alone.
PD and FD are calculated based on an ultrametric tree (hclust or phylo object).
The number and order of species in comm must be the same as in tree.
}
\examples{
comm <- matrix(c(2,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, ncol = 5, byrow = TRUE)
tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
beta.multi(comm)
beta.multi(comm, func = "Soerensen")
beta.multi(comm, tree)
beta.multi(comm, raref = 1)
beta.multi(comm, tree, T, "s", raref = 2)
}
\references{
Cardoso, P., Rigal, F., Carvalho, J.C., Fortelius, M., Borges, P.A.V., Podani, J. & Schmera, D. (2014) Partitioning taxon, phylogenetic and functional beta diversity into replacement and richness difference components. Journal of Biogeography, 41, 749-761.

Carvalho, J.C., Cardoso, P. & Gomes, P. (2012) Determining the relative roles of species replacement and species richness differences in generating beta-diversity patterns. Global Ecology and Biogeography, 21, 760-771.

Legendre, P. (2014) Interpreting the replacement and richness difference components of beta diversity. Global Ecology and Biogeography, in press.

Podani, J. & Schmera, D. (2011) A new conceptual and methodological framework for exploring and explaining pattern in presence-absence data. Oikos, 120, 1625-1638.
}

