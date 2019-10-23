
#Packages needed

load_required_packages = function(){
	library("lme4")
	library("languageR")
	library("arm")
}

# Define all input by user here, output of gene expression mega-analysis and meta-analysis functions is described at the end of function below

define_paths_and_variables = function(){

	#file location/paths
	
	#user input

	#large single gene expression data matrix created by combining all individual processed studies into one matrix
	#Column 1 name= "individualid_col" which are the unique sample IDs from all studies
	#Column 2 name= "grp" which is 2 integer groups for e.g. "1" (Group1 microarrays for e.g. "cases") and "0" (Group2 microarrays for e.g. "controls") for all samples from all studies
	#Column 3 name= "studyid_col" which is the study.name from all studies
	#Column 4 onwards = RefSeq gene symbols with residual expression values from all studies for each RefSeq gene symbols
	ExpressionFileName="/home/cluster/nodep/fseifudd/30-gene-expression-bp-data-search/09-Gene-Expression-Analysis/input/genes_all_expression_studies_transposed.txt"

	# vector - of study.names to be included in the mega-analysis (this depends on how you named your individual studies "study.name variable" in the previous individual data processing step)
	studies=c(1,3,5,7,4,19)
	
	# string - name of mega-analysis
	metaname="pfc"
	
	# integer - total number of genes that are included in the large single gene expression data matrix created by combining all individual processed studies into one matrix
	total.number.of.genes=16884

	#output

	# File = "metaname_mega_statistics.txt"
		# Columns =
		# standard error = se_fe_4dp
		# beta estimate = estimate_fe_4dp
		# fold change = foldc_fe_2dp
		# upper limit 95% = upperl_fe_2dp
		# lower limit 95% = lowerl_fe_2dp
		# pvalue = pval_fe_3sig

	# File = "metaname_meta_statistics.txt"
		# Columns =
		# standard error = se_fe_4dp
		# beta estimate = estimate_fe_4dp
		# fold change = foldc_fe_2dp
		# upper limit 95% = upperl_fe_2dp
		# lower limit 95% = lowerl_fe_2dp
		# pvalue = pval_fe_3sig

	return(list("ExpressionFileName" = ExpressionFileName, "studies" = studies, "metaname" = metaname, "total.number.of.genes" = total.number.of.genes))

}

lme4_mixed_model = function(ExpressionFileName=NULL, studies=NULL, metaname=NULL, total.number.of.genes=NULL){

	bpexpr = read.table(file=ExpressionFileName,header=T,row.names=1,colClasses=c(rep("factor",4),rep("numeric",total.number.of.genes)), na.strings="NA")
	bpexpr_sub = subset(bpexpr,studyid_col %in% studies)
	bpexpr_sub2 = subset(bpexpr_sub, select=-c(grp,studyid_col,individualid_col))
	lme4_statistics = apply(bpexpr_sub2, 2, function(x) try(lme4_mega(x, d=bpexpr_sub), silent=TRUE))

	list_lme4_statistics = lme4_statistics[sapply(lme4_statistics, function(x) !inherits(x, "try-error"))]
	table_lme4_statistics = do.call("rbind",list_lme4_statistics)
	genes = names(list_lme4_statistics)
	final_table_lme4_statistics = cbind(genes,table_lme4_statistics)
	write.table(final_table_lme4_statistics,file=paste(metaname,"_mega_statistics.txt",sep=""), sep = "\t", row.names=FALSE, quote=FALSE)
}

lm_model = function(ExpressionFileName=NULL, study=NULL, metaname=NULL, total.number.of.genes=NULL){

	bpexpr = read.table(file=ExpressionFileName,header=T,row.names=1,colClasses=c(rep("factor",4),rep("numeric",total.number.of.genes)), na.strings="NA")
	bpexpr_sub = subset(bpexpr,studyid_col==study)
	bpexpr_sub2 = subset(bpexpr_sub, select=-c(grp,studyid_col,individualid_col))
	lm_statistics = apply(bpexpr_sub2, 2, function(x) try(linear_regression(x, d=bpexpr_sub), silent=TRUE))

	list_lm_statistics = lm_statistics[sapply(lm_statistics, function(x) !inherits(x, "try-error"))]
	table_lm_statistics = do.call("rbind",list_lm_statistics)
	genes = names(list_lm_statistics)
	final_table_lm_statistics = cbind(genes,table_lm_statistics)
	write.table(final_table_lm_statistics,file=paste(metaname,"_meta_statistics.txt",sep=""), sep = "\t", row.names=FALSE, quote=FALSE)
}

# mixed model analysis for data analysis (mega-analysis)
lme4_mega = function(x,d=NULL){

	fm3 = lmer(x ~ grp + (1|studyid_col) + (1|individualid_col), data=d)
	mcmc = pvals.fnc(fm3,nsim=10000,addPlot=FALSE, ndigits=20)
	se_fe = se.fixef(fm3)[[2]]
	estimate_fe = as.numeric(mcmc$fixed[,1][2])
	pval_fe = as.numeric(mcmc$fixed[,6][2])

	foldc_fe = 2^estimate_fe
	upperl_fe = 2^(estimate_fe + (1.96 * se_fe))
	lowerl_fe = 2^(estimate_fe - (1.96 * se_fe))
	
	if(foldc_fe < 1.0){
		foldc_fe = -(1.0 / foldc_fe)
	}

	if(upperl_fe < 1.0){
		upperl_fe = -(1.0 / upperl_fe)
	}

	if(lowerl_fe < 1.0){
		lowerl_fe = -(1.0 / lowerl_fe)
	}

	se_fe_4dp = formatC(x=round(se_fe,4),digits=4,format="f")
	estimate_fe_4dp = formatC(x=round(estimate_fe,4),digits=4,format="f")
	pval_fe_3sig = format.pval(pval_fe, digits=3)
	foldc_fe_2dp = formatC(x=round(foldc_fe,2),digits=2,format="f")
	upperl_fe_2dp = formatC(x=round(upperl_fe,2),digits=2,format="f")
	lowerl_fe_2dp = formatC(x=round(lowerl_fe,2),digits=2,format="f")

	lme4_stats = cbind(se_fe_4dp, estimate_fe_4dp, foldc_fe_2dp, upperl_fe_2dp, lowerl_fe_2dp, pval_fe_3sig)
}


# linear regression analysis for data analysis (meta-analysis)
linear_regression_meta = function(x,d=NULL){
	lm_model = lm(x ~ grp, data=d)
	se_fe = coefficients(summary(lm_model))[,2][2]
	estimate_fe = coefficients(summary(lm_model))[,1][2]
	pval_fe = coefficients(summary(lm_model))[,4][2]

	foldc_fe = 2^estimate_fe
	upperl_fe = 2^(estimate_fe + (1.96 * se_fe))
	lowerl_fe = 2^(estimate_fe - (1.96 * se_fe))
	
	if(foldc_fe < 1.0){
		foldc_fe = -(1.0 / foldc_fe)
	}

	if(upperl_fe < 1.0){
		upperl_fe = -(1.0 / upperl_fe)
	}

	if(lowerl_fe < 1.0){
		lowerl_fe = -(1.0 / lowerl_fe)
	}

	se_fe_4dp = formatC(x=round(se_fe,4),digits=4,format="f")
	estimate_fe_4dp = formatC(x=round(estimate_fe,4),digits=4,format="f")
	pval_fe_3sig = format.pval(pval_fe, digits=3)
	foldc_fe_2dp = formatC(x=round(foldc_fe,2),digits=2,format="f")
	upperl_fe_2dp = formatC(x=round(upperl_fe,2),digits=2,format="f")
	lowerl_fe_2dp = formatC(x=round(lowerl_fe,2),digits=2,format="f")

	residual_stats = cbind(se_fe_4dp, estimate_fe_4dp, foldc_fe_2dp, upperl_fe_2dp, lowerl_fe_2dp, pval_fe_3sig)
}

#execute functions
load_required_packages()

variables = define_paths_and_variables()

lme4_mixed_model(ExpressionFileName = variables$ExpressionFileName, studies = variables$studies, metaname = variables$metaname, total.number.of.genes = variables$total.number.of.genes)

lme4_model(ExpressionFileName = variables$ExpressionFileName, studies = variables$studies, metaname = variables$metaname, total.number.of.genes = variables$total.number.of.genes)

