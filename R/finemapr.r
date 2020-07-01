#' Generate data for analysis in various finemapping methods
#'
#' Uses the finemapr package https://github.com/variani/finemapr
#'
#' @param region Region of the genome to extract eg 1:109317192-110317192"
#' @param id Array of GWAS studies to query. See \code{gwasinfo} for available studies
#' @param bfile If this is provided then will use the API. Default = NULL
#' @param plink_bin If null and bfile is not null then will detect packaged plink binary for specific OS. Otherwise specify path to plink binary. Default = NULL
#'
#' @export
#' @return Each id will be a list of z score data, ld matrix, and sample size
ieugwasr_to_finemapr <- function(region, id, bfile=NULL, plink_bin=NULL)
{
	id <- unique(id)
	message("Getting rsids in region")
	rsid <- ieugwasr::variants_to_rsid(region)
	message("Extracting rsids from data")
	as <- ieugwasr::associations(rsid, id, proxies=0)
	rsid_avail <- unique(as$rsid)
	message("Calculating LD for ", length(rsid_avail), " variants")
	ld <- suppressWarnings(ieugwasr::ld_matrix(rsid_avail, bfile, plink_bin, with_alleles=FALSE)) %>% greedy_remove()
	rsid_avail <- rownames(ld)
	message("Data available for ", length(rsid_avail), " variants")
	as <- subset(as, rsid %in% rsid_avail)
	out <- list()
	for(i in 1:length(unique(id)))
	{
		dat <- list()
		x <- as[as[["rsid"]] %in% rsid_avail & as[["id"]] == id[i],]
		dat[["z"]] <- dplyr::tibble(snp = x[["rsid"]], zscore = x[["beta"]] / x[["se"]])
		index <- match(x[["rsid"]], rsid_avail)
		dat[["ld"]] <- ld[index, index]
		stopifnot(all(x[["rsid"]] == rownames(dat[["ld"]])))

		n <- x[["n"]]
		if(all(is.na(n)))
		{
			g <- ieugwasr::gwasinfo(id[i])
			n <- g[["sample_size"]]
		}
		dat[["n"]] <- n
		out[[id[i]]] <- dat
	}
	class(out) <- "FinemaprList"
	return(out)
}

print.FinemaprList <- function(x)
{
	utils::str(x)
}

#' Generate data for analysis in various finemapping methods
#'
#' Uses the finemapr package https://github.com/variani/finemapr
#'
#' @param region Region of the genome to extract eg 1:109317192-110317192"
#' @param vcf list of paths to gwasvcf files
#' @param bfile If this is provided then will use the API. Default = NULL
#' @param plink_bin If null and bfile is not null then will detect packaged plink binary for specific OS. Otherwise specify path to plink binary. Default = NULL
#'
#' @importFrom purrr map
#' @importFrom gwasvcf query_gwas vcf_to_granges
#' @importFrom ieugwasr ld_matrix
#' @importFrom dplyr transmute
#' @importFrom tibble as_tibble
#' 
#' @export
#' @return Each vcf will be a list of z score data, ld matrix, and sample size
gwasvcf_to_finemapr <- function(region, vcf,
                                bfile     = NULL,
                                plink_bin = getOption("tools_plink"))
{
    out <- regions %>%
        unique() %>%
        setNames(., .) %>%
        map(~{
            ## associations
            as   <- query_gwas(vcf, chrompos = .x)
            ## filter invalids
            gr <- vcf_to_granges(as)
            rsid  <- names(gr)[ is.finite(gr$ES) &
                                is.finite(gr$SE)]
            as  <- as[rsid]
            ## get LD
            ld <- NULL
            rsid_avail <- rsid
            if(!is.null(bfile) & !is.null(plink_bin)){
                ld <- suppressWarnings({
                    ieugwasr::ld_matrix(rsid,
                                        bfile,
                                        plink_bin,
                                        with_alleles=FALSE)})
                ld <- greedy_remove(ld)
                rsid_avail <- rownames(ld)
            }
            ## subset
            as <- as[rsid_avail, ]
            ld <- ld[rsid_avail, rsid_avail]
            ## format AS
            asdf <- as %>%
                vcf_to_granges() %>%
                dplyr::as_tibble()
            ## cases
            n <-  asdf$SS
            if(all(is.na(n))){
                n <- as %>%
                    VariantAnnotation::header() %>%
                    VariantAnnotation::meta() %>%
                    .$SAMPLE %>%
                    .$TotalControls %>%
                    as.numeric()
            }
            ## put together
            dat <- list(
                z  = asdf %>%
                    transmute(rsid = ID,
                              z    = ifelse(!is.na(EZ), EZ, ES/SE)) %>%
                    as.data.frame(),
                LD       = ld,
                N        = n,
                s        = max(asdf$NC)/max(n),
                pvalues  = 10**-asdf$LP,
                MAF      = asdf$AF,
                beta     = asdf$ES,
                varbeta  = asdf$SE^2,
                type     = type,
                snp      = asdf$ID
            )
            dat
        })

    ## and return
    class(out) <- "FinemaprList"
    return(out)
}



greedy_remove <- function(ld)
{
	ind <- which(!is.finite(ld), arr.ind=TRUE)
	if(length(ind) == 0)
	{
		return(ld)
	}
	tab <- table(ind) %>% sort(decreasing=TRUE) %>% as.data.frame(stringsAsFactors=FALSE)
	rem <- c()
	for(i in 1:nrow(tab))
	{
		ind <- ind[!(ind[,1] == tab[["ind"]][i] | ind[,2] == tab[["ind"]][i]), ]
		rem <- c(rem, tab[["ind"]][i])
		if(nrow(ind) == 0) break
	}
	rem <- as.numeric(rem)
	ld <- ld[-rem, -rem]
	stopifnot(all(is.finite(ld)))
	return(ld)
}
