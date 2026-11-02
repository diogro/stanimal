#' Additive relationship matrix from a pedigree
#'
#' Computes the numerator (additive) relationship matrix \code{A} from a
#' pedigree, for use as the \code{A} argument of \code{\link{lmm_animal}}.
#'
#' @param pedigree Data frame whose first three columns are the individual, its
#'   sire and its dam, with \code{0}, \code{NA} or \code{""} for unknown
#'   parents. Individuals can be in any order, and the pedigree can contain
#'   ancestors that have no records.
#' @param Gamma Optional relationship matrix among the founders, or among
#'   founder groups when \code{groups} is given. Founders are the individuals
#'   with both parents unknown, including parents that were added because they
#'   were missing from the pedigree. The default \code{NULL} leaves founders
#'   unrelated, as in a standard pedigree. With \code{groups = NULL},
#'   \code{Gamma} must have one row and column per founder, or founder IDs as
#'   row and column names; with \code{groups}, it must have one row and column
#'   per group, or the group labels as row and column names.
#' @param groups Optional vector assigning each founder to a row and column of
#'   \code{Gamma}, for a metafounder or genetic-group model. It can be named
#'   with founder IDs, or unnamed with one entry per founder in pedigree order
#'   (individuals with both parents unknown, then any parents added as
#'   founders). When \code{Gamma} has no dimnames, groups are used in
#'   \code{sort(unique(groups))} order.
#' @return A symmetric positive-definite matrix with the individuals' IDs on both
#'   dimensions, in the order of \code{pedigree}.
#' @details The additive relationship is computed with the tabular method: with
#'   \eqn{s_i} and \eqn{d_i} the parents of individual \eqn{i},
#'   \deqn{A_{ii} = 1 + A_{s_i d_i} / 2, \qquad A_{ij} = (A_{s_i j} + A_{d_i j}) / 2,}
#'   taking a term as zero when the parent is unknown. The first line includes
#'   the inbreeding coefficient of \eqn{i} through the relationship between its
#'   parents, so pedigrees with inbreeding are handled. Individuals are
#'   processed in an order where parents come before their offspring, whatever
#'   order the pedigree is in, and a pedigree where an individual is its own
#'   ancestor is an error.
#'
#'   A parent that is not one of the individuals of the pedigree is added to
#'   \code{A} as a founder, with a warning, so a truncated pedigree or a
#'   mistyped ID does not silently drop the relationship. Subset \code{A} by
#'   name to match the rows of \code{Y}, or let \code{\link{animal_design}} do
#'   it.
#'
#'   \code{\link{animal_design}} calls this function when it is given
#'   \code{pedigree} and no genotypes.
#'
#'   With \code{Gamma}, the founder block of \code{A} is initialized with the
#'   founder relationships instead of the identity, and the same recursion
#'   propagates them. This is the \eqn{A^\Gamma} construction for related
#'   founders and metafounders (Legarra et al. 2015): \eqn{\Gamma_{ii}} is the
#'   relationship of founder \eqn{i} with itself and \eqn{\Gamma_{ij}} its
#'   relationship with founder \eqn{j}, so a founder no longer has to have
#'   \eqn{A_{ii} = 1}. When \code{groups} is given, founders in the same group
#'   share the diagonal of \code{Gamma} and different groups use its
#'   off-diagonals. A genomic relationship matrix among the founders, for
#'   example from \code{\link{genomic_relationship}}, can be used as
#'   \code{Gamma}, but it estimates realized relationships and is not the same
#'   as the base-population \eqn{\Gamma}; shrinkage, such as the \code{blend}
#'   argument of \code{genomic_relationship}, is advisable.
#'
#'   Only founders with both parents unknown are given founder relationships.
#'   An unknown parent of an individual that has one known parent is still
#'   treated as an unrelated founder, as in the standard tabular method; to
#'   model those unknown parents as part of the related base population, add
#'   them explicitly to the pedigree and assign them to a group.
#' @references Legarra, A., Christensen, O. F., Vitezica, Z. G., Aguilar, I.,
#'   & Misztal, I. (2015). Ancestral relationships using metafounders: finite
#'   ancestral populations and across population relationships. Genetics,
#'   200(2), 455-468. \doi{10.1534/genetics.115.177014}
#' @seealso \code{\link{genomic_relationship}}, \code{\link{animal_design}}
#' @examples
#' ped <- data.frame(id = 1:6, sire = c(0, 0, 1, 1, 4, 4),
#'                   dam = c(0, 0, 2, 2, 3, 5))
#' A <- makeA(ped)
#' round(A, 3)
#' # full-sib mating: individual 6 is the offspring of the sibs 4 and 5
#' ped2 <- data.frame(id = 1:6, sire = c(NA, NA, NA, 1, 1, 4),
#'                    dam = c(NA, NA, NA, 2, 2, 5))
#' diag(makeA(ped2))   # 1.25 for the inbred individual
#'
#' # related founders: founders 1 and 2 have relationship 0.5
#' Gamma <- matrix(c(1, 0.5, 0.5, 1), 2,
#'                 dimnames = list(c("1", "2"), c("1", "2")))
#' ped3 <- data.frame(id = 1:4, sire = c(NA, NA, 1, 1),
#'                    dam = c(NA, NA, 2, 2))
#' A3 <- makeA(ped3, Gamma = Gamma)
#' round(A3, 3)
#'
#' # founder groups: 1 and 2 in group A, 3 in group B
#' Gamma_groups <- matrix(c(1, 0.2, 0.2, 1), 2,
#'                        dimnames = list(c("A", "B"), c("A", "B")))
#' ped4 <- data.frame(id = 1:5, sire = c(NA, NA, NA, 1, 1),
#'                    dam = c(NA, NA, NA, 2, 3))
#' round(makeA(ped4, Gamma = Gamma_groups, groups = c("A", "A", "B")), 3)
#' @export
makeA <- function(pedigree, Gamma = NULL, groups = NULL) {
  ped <- .clean_pedigree(pedigree)
  ids <- as.character(ped$id)
  if (anyDuplicated(ids))
    stop("`pedigree` has duplicated individuals.", call. = FALSE)

  parent_ids <- unique(c(as.character(ped$sire), as.character(ped$dam)))
  parent_ids <- parent_ids[!is.na(parent_ids)]
  added <- setdiff(parent_ids, ids)
  if (length(added))
    warning("parents that are not individuals in the pedigree were added as ",
            "founders: ", paste(utils::head(added, 5L), collapse = ", "),
            if (length(added) > 5L) ", ..." else "", call. = FALSE)

  all_ids <- c(ids, added)
  n <- length(all_ids)
  # founders added for missing parents have no parents of their own
  founder <- rep(0L, length(added))
  sire <- c(.parent_positions(ped$sire, all_ids), founder)
  dam <- c(.parent_positions(ped$dam, all_ids), founder)
  order <- .pedigree_order(sire, dam, n)

  founders <- which(sire == 0L & dam == 0L)
  A <- diag(n)
  if (!is.null(Gamma)) {
    A[founders, founders] <- .founder_relationships(
      Gamma, groups, all_ids[founders])
  } else if (!is.null(groups)) {
    stop("`groups` requires `Gamma`.", call. = FALSE)
  }

  processed <- integer(n)
  count <- 0L
  for (position in seq_len(n)) {
    i <- order[position]
    # Founders already carry Gamma in their block; filling their row would
    # overwrite those relationships with zeros.
    if (count && (sire[i] > 0L || dam[i] > 0L)) {
      previous <- processed[seq_len(count)]
      from_sire <- if (sire[i] > 0L) A[sire[i], previous] else 0
      from_dam <- if (dam[i] > 0L) A[dam[i], previous] else 0
      values <- 0.5 * (from_sire + from_dam)
      A[i, previous] <- values
      A[previous, i] <- values
    }
    if (sire[i] > 0L && dam[i] > 0L)
      A[i, i] <- 1 + 0.5 * A[sire[i], dam[i]]
    count <- count + 1L
    processed[count] <- i
  }

  dimnames(A) <- list(all_ids, all_ids)
  A
}

# The first three columns of a pedigree, named id, sire and dam, with unknown
# parents coded as NA.
.clean_pedigree <- function(pedigree) {
  ped <- as.data.frame(pedigree, stringsAsFactors = FALSE)
  if (ncol(ped) < 3L)
    stop("`pedigree` must have at least three columns: individual, sire and dam.",
         call. = FALSE)
  ped <- ped[, 1:3, drop = FALSE]
  names(ped) <- c("id", "sire", "dam")
  ped[] <- lapply(ped, function(x) if (is.factor(x)) as.character(x) else x)
  if (anyNA(ped$id))
    stop("`pedigree` has individuals with a missing ID.", call. = FALSE)
  for (parent in c("sire", "dam")) {
    x <- ped[[parent]]
    unknown <- if (is.numeric(x)) !is.na(x) & x == 0
               else !is.na(x) & x %in% c("", "0")
    x[unknown] <- NA
    ped[[parent]] <- x
  }
  ped
}

# Positions of the parents in the ID vector, with 0 for an unknown parent.
.parent_positions <- function(x, ids) {
  position <- match(as.character(x), ids)
  position[is.na(position)] <- 0L
  position
}

# Expands a founder-by-founder or group-by-group relationship matrix to the
# founders of the pedigree, in the order of `founders`.
.founder_relationships <- function(Gamma, groups, founders) {
  if (!is.matrix(Gamma) || !is.numeric(Gamma) || any(!is.finite(Gamma)))
    stop("`Gamma` must be a finite numeric matrix.", call. = FALSE)
  if (nrow(Gamma) != ncol(Gamma))
    stop("`Gamma` must be a square matrix.", call. = FALSE)
  if (nrow(Gamma) < 1L)
    stop("`Gamma` must have at least one row and column.", call. = FALSE)
  Gamma <- .check_symmetric(Gamma, "Gamma")
  values <- .relationship_eigenvalues(Gamma)
  if (min(values) <= .eigen_tolerance(values))
    stop("`Gamma` must be positive definite.", call. = FALSE)

  gamma_ids <- rownames(Gamma)
  if (is.null(gamma_ids)) gamma_ids <- colnames(Gamma)
  if (!is.null(gamma_ids)) {
    if (anyNA(gamma_ids) || anyDuplicated(gamma_ids))
      stop("`Gamma` row and column names must be unique and not missing.",
           call. = FALSE)
    if (!is.null(rownames(Gamma)) && !is.null(colnames(Gamma)) &&
        !identical(rownames(Gamma), colnames(Gamma)))
      stop("`Gamma` must have the same row and column names in the same order.",
           call. = FALSE)
    dimnames(Gamma) <- list(gamma_ids, gamma_ids)
  }

  if (is.null(groups)) {
    if (!is.null(gamma_ids)) {
      missing <- setdiff(founders, gamma_ids)
      if (length(missing))
        stop("`Gamma` has no row for founder(s): ",
             paste(utils::head(missing, 5L), collapse = ", "), call. = FALSE)
      return(Gamma[founders, founders, drop = FALSE])
    }
    if (nrow(Gamma) != length(founders))
      stop("`Gamma` must have one row and column per founder (",
           length(founders), ") or founder IDs as dimnames.", call. = FALSE)
    return(Gamma)
  }

  if (!is.null(names(groups))) {
    if (anyDuplicated(names(groups)))
      stop("`groups` names must be unique.", call. = FALSE)
    missing <- setdiff(founders, names(groups))
    if (length(missing))
      stop("`groups` has no entry for founder(s): ",
           paste(utils::head(missing, 5L), collapse = ", "), call. = FALSE)
    groups <- groups[founders]
  } else if (length(groups) != length(founders)) {
    stop("`groups` must have one entry per founder (", length(founders),
         ") or be named with founder IDs.", call. = FALSE)
  }
  if (is.factor(groups)) {
    labels <- levels(droplevels(groups))
    groups <- as.character(groups)
  } else {
    labels <- as.character(sort(unique(groups)))
    groups <- as.character(groups)
  }
  if (anyNA(groups))
    stop("`groups` must not contain missing values.", call. = FALSE)

  if (!is.null(gamma_ids)) {
    missing <- setdiff(labels, gamma_ids)
    if (length(missing))
      stop("`Gamma` has no row for group(s): ",
           paste(utils::head(missing, 5L), collapse = ", "), call. = FALSE)
    Gamma <- Gamma[labels, labels, drop = FALSE]
  } else if (nrow(Gamma) != length(labels)) {
    stop("`Gamma` must have one row and column per group (", length(labels),
         ") or group names as dimnames.", call. = FALSE)
  }
  Gamma[match(groups, labels), match(groups, labels), drop = FALSE]
}

# Order in which individuals can be processed so that parents come before their
# offspring.
.pedigree_order <- function(sire, dam, n) {
  n_parents <- (sire > 0L) + (dam > 0L)
  children <- vector("list", n)
  for (i in seq_len(n)) {
    for (parent in c(sire[i], dam[i])) {
      if (parent > 0L) children[[parent]] <- c(children[[parent]], i)
    }
  }

  ready <- which(n_parents == 0L)
  order <- integer(n)
  filled <- 0L
  head <- 1L
  while (head <= length(ready)) {
    i <- ready[head]
    head <- head + 1L
    filled <- filled + 1L
    order[filled] <- i
    for (child in children[[i]]) {
      n_parents[child] <- n_parents[child] - 1L
      if (n_parents[child] == 0L) ready <- c(ready, child)
    }
  }
  if (filled < n)
    stop("`pedigree` has a loop: some individuals appear as their own ",
         "ancestors.", call. = FALSE)
  order
}
