library(phyloseq)
library(vegan)

adonis_table <- function(phy, d, ...){
	subs <- function(dist, index) as.dist(as.matrix(dist)[index,index])
	contcheck <- function(var){
		suppressWarnings(
		if(all(is.na(as.numeric(var))) | length(unique(var))<10){
			return(var)
		} else {
			return(as.numeric(var))
		}
		)
	}
	suppressWarnings(
		tab <- phy %>%
			get_variable(.) %>%
			gather(variable, value, ...) 
		)
	tab %>% 
		group_by(variable) %>%
		do( adonis(subs(d,!is.na(.$value)) ~ contcheck(value), .[!is.na(.$value),]) %>% 
			( function(x) {data.frame(	x$aov.tab[1,4:6], 
										df = x$aov.tab[1,1],
										n = nrow(x$model.matrix),
										perms = nrow(x$f.perms)
										)} )
			) 	 %>%
		ungroup %>%
		rename(pval = `Pr..F.`)
}

# Examples with phyloseq dummy data
data(soilrep)
adonis_table(soilrep, distance(soilrep, "bray"), Treatment, warmed, clipped, Sample)
data(GlobalPatterns)
adonis_table(GlobalPatterns, distance(GlobalPatterns, "bray"), SampleType, Description)
