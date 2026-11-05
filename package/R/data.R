#' Morphological and life-history data for five mouse lines
#'
#' Cranial measurements and individual information for mice from the long-term
#' artificial selection experiment on body weight analyzed by Penna et al.
#' (2017). Selection for body weight at 49 days of age was carried out
#' downwards in the \code{s} and \code{h} lines and upwards in the \code{sp} and
#' \code{hp} lines, alongside a randomly mated control line \code{t}; the
#' \code{s}/\code{sp} lines were inbred by limiting population size and the
#' \code{h}/\code{hp} lines by full-sib mating in the first generations.
#'
#' Every individual carries a full set of 35 cranial Euclidean distances,
#' measured as the mean of two replicas, plus the catalogue numbers of its
#' parents. The
#' distances are moderately integrated: the first principal component of their
#' correlation matrix explains about 40\% of the variance and the mean pairwise
#' correlation is 0.35, so the traits carry partly independent information and
#' a handful of them is enough for a well-identified model.
#'
#' @format A data frame with 329 rows and 46 columns. Row names are the
#'   individual IDs (\code{as.character(ID)}).
#'   \describe{
#'     \item{\code{.id}}{line label combining the selection direction and the
#'       line, e.g. \code{"control.t"};}
#'     \item{\code{ID}}{individual catalogue number;}
#'     \item{\code{SEX}}{sex, \code{"F"} or \code{"M"};}
#'     \item{\code{GER}}{generation of the long-term experiment;}
#'     \item{\code{LIN}}{experimental line: \code{t} (control), \code{s} and
#'       \code{h} (downward selection), \code{sp} and \code{hp} (upward
#'       selection);}
#'     \item{\code{MADRE}}{catalogue number of the dam;}
#'     \item{\code{PADRE}}{catalogue number of the sire;}
#'     \item{\code{AGE}}{age at sacrifice, in days;}
#'     \item{\code{P49}}{body weight at 49 days of age;}
#'     \item{\code{line}}{the line as a factor with levels \code{t},
#'       \code{h'}, \code{s'}, \code{h} and \code{s};}
#'     \item{\code{selection}}{direction of selection: \code{"control"},
#'       \code{"downwards"} or \code{"upwards"};}
#'     \item{\code{IS_PM}, \code{IS_NSL}, \code{IS_PNS}, \code{PM_ZS},
#'       \code{PM_ZI}, \code{PM_MT}, \code{NSL_NA}, \code{NSL_ZS},
#'       \code{NSL_ZI}, \code{NA_BR}, \code{NA_PNS}, \code{BR_PT},
#'       \code{BR_APET}, \code{PT_APET}, \code{PT_BA}, \code{PT_EAM},
#'       \code{PT_ZYGO}, \code{PT_TSP}, \code{ZS_ZI}, \code{ZI_MT},
#'       \code{ZI_ZYGO}, \code{ZI_TSP}, \code{MT_PNS}, \code{PNS_APET},
#'       \code{APET_BA}, \code{APET_TS}, \code{BA_EAM}, \code{EAM_ZYGO},
#'       \code{ZYGO_TSP}, \code{LD_AS}, \code{BR_LD}, \code{OPI_LD},
#'       \code{PT_AS}, \code{JP_AS}, \code{BA_OPI}}{35 cranial Euclidean
#'       distances, named after the landmarks they span. They are on the
#'       measured scale; \code{log()} is the usual transformation before
#'       fitting.}
#'   }
#' @source Penna et al. (2017), archived at Dryad under CC0 1.0
#'   (\doi{10.5061/dryad.5gr8r}) and
#'   distributed by the MIT-licensed companion
#'   package \url{https://github.com/diogro/ratones}.
#' @references Penna, A., Melo, D., Bernardi, S., Oyarzabal, M. I., & Marroig,
#'   G. (2017). The evolution of phenotypic integration: How directional
#'   selection reshapes covariation in mice. Evolution, 71(10), 2370-2380.
#'   \doi{10.1111/evo.13304}
#' @section Weak identification:
#'   This is a long-term artificial-selection experiment that happens to carry a
#'   pedigree, not a design built to estimate a genetic covariance matrix. The
#'   animal model fits, but 329 individuals from 58 sires identify a 35-trait
#'   \eqn{\mathbf{G}} only weakly: the split between the genetic and the
#'   environmental covariance is left largely to the priors, so a heritability
#'   read from \code{\link{lmm_animal}} on these traits depends on the prior as
#'   much as on the data, and can move by a factor of eight between defensible
#'   choices. Use it as the small worked example it is, and use
#'   \code{\link{opossum}}, which carries marker genotypes and a design that does
#'   identify \eqn{\mathbf{G}}, for anything about genetic covariance.
#' @seealso \code{\link{ratones_pedigree}}, \code{\link{animal_design}},
#'   \code{\link{lm_phen_cov}}
#' @examples
#' data(ratones_morpho)
#' dim(ratones_morpho)
#' table(ratones_morpho$selection)
#' traits <- names(ratones_morpho)[12:46]
#' pc1_variance_share(log(as.matrix(ratones_morpho[, traits])))
"ratones_morpho"

#' Pedigree of the mouse lines
#'
#' The pedigree behind \code{\link{ratones_morpho}}: all 525 individuals of the
#' experiment, including the ancestors that were never measured. Used to build
#' the additive relationship matrix \code{A} for \code{\link{lmm_animal}}.
#'
#' The pedigree is complete and consistent, and the animal model will run on it,
#' but it does not make the mice a quantitative-genetic design. Read the weak
#' identification section of \code{\link{ratones_morpho}} before using it for
#' anything beyond a worked example, and prefer \code{\link{opossum}} and its
#' genomic relationship matrix \code{opossum_grm} when the genetic covariance is
#' the object of inference.
#'
#' @format A data frame with 525 rows and 3 columns, in the order expected by
#'   \code{\link{animal_design}}:
#'   \describe{
#'     \item{\code{id}}{individual catalogue number;}
#'     \item{\code{sire}}{catalogue number of the sire, \code{NA} for a
#'       founder or an unknown parent;}
#'     \item{\code{dam}}{catalogue number of the dam, \code{NA} for a founder
#'       or an unknown parent.}
#'   }
#' @source Penna et al. (2017), archived at Dryad under CC0 1.0
#'   (\doi{10.5061/dryad.5gr8r}) and
#'   distributed by the MIT-licensed companion
#'   package \url{https://github.com/diogro/ratones}.
#' @references Penna, A., Melo, D., Bernardi, S., Oyarzabal, M. I., & Marroig,
#'   G. (2017). The evolution of phenotypic integration: How directional
#'   selection reshapes covariation in mice. Evolution, 71(10), 2370-2380.
#'   \doi{10.1111/evo.13304}
#' @seealso \code{\link{ratones_morpho}}, \code{\link{makeA}}
#' @examples
#' data(ratones_pedigree)
#' dim(ratones_pedigree)
#' head(ratones_pedigree)
#' sum(is.na(ratones_pedigree$sire))  # founders
"ratones_pedigree"

#' Opossum cranial measurements and a genomic relationship matrix
#'
#' Cranial Euclidean distances and individual information for opossums from the
#' colony of the South Texas Diabetes and Obesity Institute, analysed by
#' \cite{Porto et al. (2016)}. The traits are the same 35 landmarks as
#' \code{\link{ratones_morpho}}, so the two datasets can be compared directly.
#'
#' Unlike the mice, these animals carry marker genotypes, so the relationship
#' matrix is genomic rather than genealogical. That matters for these data: the
#' colony pedigree shipped with the original archive gives every pair a kinship
#' coefficient between 0.69 and 0.93, which leaves additive variance
#' indistinguishable from a shared colony effect. The genomic matrix separates
#' relatives from non-relatives and makes the genetic covariance identifiable.
#'
#' The distances are stored as measured. \cite{Porto et al. (2016)} analysed them
#' untransformed, and the phenomenon is confirmed by their published phenotypic
#' covariance, whose diagonal correlates 0.834 with the variances of the raw
#' distances against 0.046 with those of their logarithms.
#'
#' @format Two objects.
#'   \describe{
#'     \item{\code{opossum}}{a data frame with 553 rows and 38 columns:
#'       \code{id}, \code{sex}, \code{age} in days, and the 35 cranial distances
#'       as measured, named after the landmarks they span, with underscores as in
#'       \code{\link{ratones_morpho}}. Row names are the individual IDs. 22
#'       individuals are missing every trait, which \code{\link{lmm_animal}}
#'       and \code{\link{lm_phen_cov}} handle, and 505 have no missing value in
#'       any trait or covariate.}
#'     \item{\code{opossum_grm}}{a 553 by 553 genomic relationship matrix
#'       computed as \eqn{ZZ^\top/m} from 3,696 markers, \eqn{Z} being the
#'       marker matrix centred by twice the allele frequency and \eqn{m} the
#'       number of markers. Rows and columns are named with the individual IDs
#'       of \code{opossum}, so the rows of the subset being analysed can be
#'       passed as the \code{A} argument of \code{\link{lmm_animal}} directly,
#'       as in the example. The full matrix is only positive semidefinite, so
#'       pass \code{blend_A = TRUE} if it is given whole. The diagonal runs from
#'       0.22 to 0.59, as expected for a centred marker matrix rather than a
#'       pedigree.}
#'   }
#' @source Porto, A., Schmelter, R., VandeBerg, J. L., Marroig, G., & Cheverud,
#'   J. M. (2016). Evolution of the genotype-to-phenotype map and the cost of
#'   pleiotropy in mammals. Genetics, 204(4), 1601-1612.
#'   \doi{10.1534/genetics.116.189431}. Data deposited on figshare as
#'   \dQuote{Data - Genetics - Porto et al 2016},
#'   \url{https://figshare.com/articles/dataset/4055961}.
#'
#'   The coded marker files in that archive carry no individual identifiers, so
#'   the animals are matched to the cranial measurements by row order. The two
#'   sources agree in length at 553, and the identities are attached here as the
#'   row and column names of \code{opossum_grm}.
#' @references Porto, A., Schmelter, R., VandeBerg, J. L., Marroig, G., &
#'   Cheverud, J. M. (2016). Evolution of the genotype-to-phenotype map and the
#'   cost of pleiotropy in mammals. Genetics, 204(4), 1601-1612.
#'   \doi{10.1534/genetics.116.189431}
#' @seealso \code{\link{ratones_morpho}}, \code{\link{lmm_animal}}
#' @aliases opossum_grm
#' @examples
#' data(opossum)
#' dim(opossum)
#' dim(opossum_grm)
#' traits <- setdiff(names(opossum), c("id", "sex", "age"))
#' \dontrun{
#' # Phenotypic covariance of four traits, no pedigree needed
#' lm_phen_cov(as.matrix(opossum[, traits[1:4]]),
#'             model.matrix(~ sex + age, opossum))
#'
#' # Genomic animal model on the same four traits
#' design <- animal_design(
#'   stats::reformulate(c("sex", "age"),
#'                      response = sprintf("cbind(%s)", paste(traits[1:4],
#'                                                            collapse = ", "))),
#'   data = opossum, pedigree = NULL)
#' lmm_animal(design$Y, design$X, opossum_grm[rownames(design$Y),
#'                                            rownames(design$Y)])
#' }
"opossum"
