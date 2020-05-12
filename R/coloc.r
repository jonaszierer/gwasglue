#' Generate coloc dataset from vcf files
#'
#' @param vcf1 VCF object or path to vcf file
#' @param vcf2 VCF object or path to vcf file
#' @param chrompos Character of chr:pos1-pos2
#'
#' @export
#' @return List of datasets to feed into coloc
gwasvcf_to_coloc <- function(vcf1, vcf2, chrompos)
{
	## TODO: binary or quantitative traits
	## TODO: multiallelic variants

	o <- gwasvcf::vcflist_overlaps(list(vcf1, vcf2), chrompos)
	vcf1 <- o[[1]]
	vcf2 <- o[[2]]

	if(length(vcf1) == 0)
	{
		message("No overlaps in vcf1")
		return(NULL)
	}
	if(length(vcf2) == 0)
	{
		message("No overlaps in vcf2")
		return(NULL)
	}

	stopifnot(length(vcf1) == length(vcf2))
	tab1 <- vcf1 %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble()
	tab2 <- vcf2 %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble()
	index <- as.character(tab1$REF) == as.character(tab2$REF) &
			as.character(tab1$ALT) == as.character(tab2$ALT) &
			as.character(tab1$seqnames) == as.character(tab2$seqnames) &
			tab1$start == tab2$start
	stopifnot(sum(index) > 0)

	type1 <- ifelse(VariantAnnotation::header(vcf1) %>%
		VariantAnnotation::meta() %>%
		{.[["SAMPLE"]][["StudyType"]]} == "Continuous", "quant", "cc")

	type2 <- ifelse(VariantAnnotation::header(vcf2) %>%
		VariantAnnotation::meta() %>%
		{.[["SAMPLE"]][["StudyType"]]} == "Continuous", "quant", "cc")

	out1 <- tab1[index,] %>% {list(pvalues = 10^-.$LP, N = .$SS, MAF = .$AF, beta = .$ES, varbeta = .$SE^2, type = type1, snp = names(vcf2)[index], z = .$ES / .$SE, chr = .$seqnames, pos = .$start, id = VariantAnnotation::samples(VariantAnnotation::header(vcf1))[1])}
	out2 <- tab2[index,] %>% {list(pvalues = 10^-.$LP, N = .$SS, MAF = .$AF, beta = .$ES, varbeta = .$SE^2, type = type2, snp = names(vcf2)[index], z = .$ES / .$SE, chr = .$seqnames, pos = .$start, id = VariantAnnotation::samples(VariantAnnotation::header(vcf2))[1])}

	if(type1 == "cc")
	{
		out1$s <- mean(tab1$NC / tab1$SS, na.rm=TRUE)
	} else 

	if(type2 == "cc")
	{
		out2$s <- mean(tab1$NC / tab1$SS, na.rm=TRUE)
	}

	return(list(dataset1=out1, dataset2=out2))
}




#' Generate coloc v4  dataset from vcf files
#'
#' @param vcf1 VCF object or path to vcf file
#' @param vcf2 VCF object or path to vcf file
#' @param chrompos Character of chr:pos1-pos2
#'
#' @export
#' @return List of datasets to feed into coloc (v4)
gwasvcf_to_coloc4 <- function(vcf1, vcf2, chrompos,
                              bfile,
                              plink_bin = getOption("tools_plink")){
    
    data <- list(d1 = vcf1, d2 = vcf2) %>% map(~{
        ## associations
        query_gwas(.x, chrompos = chrompos)
    })
    
    ## exclude SNPs with mismatchd alleles
    alleles  <- data %>%
        map(~{
            .x %>%
                vcf_to_granges() %>%
                dplyr::as_tibble() %>%
                transmute(ID = names(.x), start ,end, REF, ALT)
        })
    alleles <- full_join(alleles[[1]],
                         alleles[[2]],
                         by = "ID", suffix = c(".1", ".2")) 
    excl_alleles <- alleles %>%
        mutate(excl = #(start.1 != start.2) |
                   #(end.1 != end.2) |
                   ((REF.1 != REF.2) & (REF.1 != ALT.2)) |
                   ((ALT.1 != ALT.2) & (ALT.1 != REF.2))) %>% 
        filter(!is.na(excl) & excl) %>%
        pull(ID)
    
    ## LD
    rs <- data %>%
        map(names) %>%
        unlist() %>%
        unique()
    ld <- suppressWarnings({
        ieugwasr::ld_matrix(rs,
                            bfile,
                            plink_bin,
                            with_alleles=FALSE)})

    ## exclude SNPs missing from LD
    excl_ld <- rs[ !(rs %in% rownames(ld)) ]
    
    ## FORMAT
    data <- data %>% map(~{
        ## filter invalids
        gr   <- vcf_to_granges(.x)
        rsid <- names(gr)[ is.finite(gr$ES) &
                            is.finite(gr$SE)]
        rsid <- rsid[ !(rsid %in% excl_alleles)]
        rsid <- rsid[ !(rsid %in% excl_ld)]
        as   <- .x[rsid]
        ## type
        type <- as %>%
            VariantAnnotation::header() %>%
            VariantAnnotation::meta() %>%
            .$SAMPLE %>%
            .$StudyType
        type <- ifelse(type == "Continuous", "quant", "cc")
        ## data.frame
        asdf <- as %>%
            vcf_to_granges() %>%
            dplyr::as_tibble() %>%
            mutate(ID = names(as))
        ## MAF
        if(all(is.na(asdf$AF))){
            MAF <- af_matrix(variants  = asdf$ID,
                             bfile     = bfile,
                             plink_bin = plink_bin) %>%
                as.data.frame() %>%
                dplyr::rename(AF = MAF) %>%
                rownames_to_column("ID")
            asdf <- asdf %>%
                dplyr::select(-AF) %>%
                left_join(MAF, by = "ID")
        }
        ## put together
        list(N        = asdf$SS,
             s        = asdf$NC/asdf$SS,
             pvalues  = 10**-asdf$LP,
             MAF      = asdf$AF,
             beta     = asdf$ES,
             varbeta  = asdf$SE^2,
             type     = type,
             position = asdf$start,
             snp      = asdf$ID)
    })

    ## RETURN
    data$LD <- ld
    return(data)
}



#' Generate coloc dataset from the IEU GWAS database
#'
#' @param id1 ID for trait 1
#' @param id2 ID for trait 2
#' @param chrompos Character of chr:pos1-pos2
#' @param type1 Provide "cc" or "quant" to override automatic lookup of trait type for trait 1
#' @param type2 Provide "cc" or "quant" to override automatic lookup of trait type for trait 2
#'
#' @export
#' @return List of datasets to feed into coloc
ieugwasr_to_coloc <- function(id1, id2, chrompos, type1=NULL, type2=NULL)
{
	tab1 <- ieugwasr::associations(id=id1, variants=chrompos) %>% subset(., !duplicated(rsid))
	tab2 <- ieugwasr::associations(id=id2, variants=chrompos) %>% subset(., !duplicated(rsid))
	commonsnps <- tab1$rsid[tab1$rsid %in% tab2$rsid]
	tab1 <- tab1[tab1$rsid %in% commonsnps, ] %>% arrange(rsid)
	tab2 <- tab2[tab2$rsid %in% commonsnps, ] %>% arrange(rsid)
	stopifnot(all(tab1$rsid == tab2$rsid))

	index <- as.character(tab1$ea) == as.character(tab2$ea) &
			as.character(tab1$nea) == as.character(tab2$nea) &
			as.character(tab1$rsid) == as.character(tab2$rsid) &
			tab1$position == tab2$position
	stopifnot(sum(index) > 0)
	tab1$eaf <- as.numeric(tab1$eaf)
	tab2$eaf <- as.numeric(tab2$eaf)
	tab1$eaf[which(tab1$eaf > 0.5)] <- 1 - tab1$eaf[which(tab1$eaf > 0.5)]
	tab2$eaf[which(tab2$eaf > 0.5)] <- 1 - tab2$eaf[which(tab2$eaf > 0.5)]

	info1 <- ieugwasr::gwasinfo(id1)
	type1 <- get_type(info1, type1)
	info2 <- ieugwasr::gwasinfo(id2)
	type2 <- get_type(info2, type2)

	tab1 <- tab1[index,] %>% {list(pvalues = .$p, N = .$n, MAF = .$eaf, beta = .$beta, varbeta = .$se^2, type = type1, snp = .$rsid, z = .$beta / .$se, chr = .$chr, pos = .$position, id = id1)}
	tab2 <- tab2[index,] %>% {list(pvalues = .$p, N = .$n, MAF = .$eaf, beta = .$beta, varbeta = .$se^2, type = type2, snp = .$rsid, z = .$beta / .$se, chr = .$chr, pos = .$position, id = id2)}

	if(type1 == "cc")
	{
		tab1$s <- info1$ncase / info1$sample_size
	} else 

	if(type2 == "cc")
	{
		tab2$s <- info2$ncase / info2$sample_size
	}

	return(list(dataset1=tab1, dataset2=tab2))
}


get_type <- function(info, typex)
{
	if(!is.null(typex))
	{
		stopifnot(! typex %in% c("cc", "quant"))
		return(typex)
	} else if(is.na(info$unit)) {
		message("Type information not available for ", info$id, ". Assuming 'quant' but override using 'type' arguments.")
		return("quant")
	} else {
		return(ifelse(.$unit %in% c("logOR", "log odds"), "cc", "quant"))
	}
}

