
# Packages needed

load_required_packages = function(){
	library(affy)
	library(frma)
	library(panp)
	library(sva)
	library(lme4)
	library(languageR)
	library(arm)
}


# Define all input by user here, output of gene expression data processing script is described at the end of function below

define_paths_and_variables = function(){

	#file location/paths
	
	#user input

	#if providing CEL files define directory path to CEL files here and define the normalized.gedm.file.path variable below to normalized.gedm.file.path = NULL
	CELfile.directory.path = NULL

	#normalized gene expression data matrix (gedm) (rows = probes/probe sets/genes/features, columns = samples)
	#header needed
	normalized.gedm.file.path = "/home/cluster/nodep/fseifudd/30-gene-expression-bp-data-search/09-Gene-Expression-Analysis/input/study1_frma_ac.txt"
	
	#Column1 with Unique Identifier (ID) matching the column header in your gene expression data matrix (gedm) or CEL file names if providing CEL files
	#Column2 with 2 integer groups for e.g. "1" (Group1 microarrays for e.g. "cases") and "0" (Group2 microarrays for e.g. "controls")
	#header needed
	DemographicFileName.path = "/home/cluster/nodep/fseifudd/30-gene-expression-bp-data-search/09-Gene-Expression-Analysis/input/demoAnnTable1_bp_controls.txt"
	
	#download additional input files provided with script for hgu133a, hgu133plus2 and hgu95av2 microarray platforms and provide path below
	#Using Affymetrix annotation tables, we identified all probesets labeled as Negative Strand Matching Probesets that were not characterized with “cross hyb,” indicating the probeset may match to another gene.
	#We then mapped the NSMPs to the ENSEMBL transcript database version GRCh37.p8 and classified the NSMPs into four categories;
	#(i) probesets that did not map to a transcript at all (used for presence/absence calling);
	#(ii) probesets that detected sense transcripts;
	#(iii) probesets that detected antisense transcripts and
	#(iv) probesets that detected a sense transcript that overlaps with an antisense transcript
	NSMP.no.transcript.file.path = "/home/cluster/nodep/fseifudd/30-gene-expression-bp-data-search/09-Gene-Expression-Analysis/hgu133/HG-U133A.na32.annot.Negative.Strand.Matching.Probes.final.csv_merged_with_ENSEMBLE_alignment_to_transcripts_strand_info_nsmps_no_transcripts.txt"
	NSMP.sense.transcript.file.path = "/home/cluster/nodep/fseifudd/30-gene-expression-bp-data-search/09-Gene-Expression-Analysis/hgu133/HG-U133A.na32.annot.Negative.Strand.Matching.Probes.final.csv_merged_with_ENSEMBLE_alignment_to_transcripts_strand_info_sense_transcripts.txt"
	NSMP.antisense.transcript.file.path = "/home/cluster/nodep/fseifudd/30-gene-expression-bp-data-search/09-Gene-Expression-Analysis/hgu133/HG-U133A.na32.annot.Negative.Strand.Matching.Probes.final.csv_merged_with_ENSEMBLE_alignment_to_transcripts_strand_info_antisense_transcripts.txt"
	NSMP.overlap.transcript.file.path = "/home/cluster/nodep/fseifudd/30-gene-expression-bp-data-search/09-Gene-Expression-Analysis/hgu133/HG-U133A.na32.annot.Negative.Strand.Matching.Probes.final.csv_merged_with_ENSEMBLE_alignment_to_transcripts_strand_info_overlap_transcripts.txt"
	
	#download Presence/Absence calling script provided by
	#Warren P: panp: Presence-Absence Calls from Negative Strand Matching Probesets. R package version 1.26.0. 2007
	PANP.script.path = "/home/cluster/nodep/fseifudd/30-gene-expression-bp-data-search/09-Gene-Expression-Analysis/panp.R"
	
	#download additional Jetset probeset to gene mapping file provided with script for hgu133a, hgu133plus2 and hgu95av2 microarray platforms by -
	#Li Q, Birkbak NJ, Gyorffy B, Szallasi Z, Eklund AC: Jetset: selecting the optimal microarray probe set to represent a gene. BMC Bioinformatics 2011, 12:474-2105-12-474.
	#Column1 with probeset_id
	#Column2 with gene_symbol (for e.g. from RefSeq)
	#header needed
	probeset.to.gene.mapping.file.path = "/home/cluster/nodep/fseifudd/30-gene-expression-bp-data-search/09-Gene-Expression-Analysis/jetset/jetset.scores.hgu133a_1.2.0.best.probeset.to.gene.txt"

	#variables

	# vector of strings for outliers to be excluded from the gedm, comment line below if ther are no outliers
	outlier.samples = c("A-19","A-28","A-29","A-31","A-41","A-50","A-59","A-87")

	# strings

	# study name
	study.name = "AltarA"

	# microarray platorm being analyzed for e.g "hgu133a", "hgu133plus2", "hgu95av2", etc...
	arrayplatform = "hgu133a"
	
	return(list("CELfile.directory.path" = CELfile.directory.path, "normalized.gedm.file.path"=normalized.gedm.file.path, "DemographicFileName.path"=DemographicFileName.path, "outlier.samples" =outlier.samples, "study.name"=study.name, "arrayplatform"=arrayplatform, "NSMP.no.transcript.file.path"=NSMP.no.transcript.file.path, "NSMP.sense.transcript.file.path"=NSMP.sense.transcript.file.path, "NSMP.antisense.transcript.file.path"=NSMP.antisense.transcript.file.path, "NSMP.overlap.transcript.file.path"=NSMP.overlap.transcript.file.path, "PANP.script.path"=PANP.script.path, "probeset.to.gene.mapping.file.path"=probeset.to.gene.mapping.file.path))

	#output
	
	###if your input is CEL files you will get "additional output" as an Affy RObject (for later use) and a normalized gene expression data matrix
	# RObject = "study.name_affy.RObject.R"
	# File = "study.name_fRMA.normalized.txt"
	#	# normalized gene expression data matrix (gedm) (rows = probes/probe sets/genes/features, columns = samples)
	# Plot = "study.name_boxplots.pre.and.post.normalization.png"
	#	# box plots showing distribution of intensities for each microarrays/samples pre and post normalization
	###

	# Plot = "study.name_boxplots.post.normalization.png"
		# box plots showing distribution of intensities for each microarrays/samples post normalization
	# Plot = "study.name_NegativeStrandMatchingProbesets.intensity.distributions.png"
		# density plot showing distribution of intensities for different categories of Negative Strand Matching Probesets (PANP provided NSMPs, all probes, no.transcript, sense.transcript, antisense.transcript, overlap.transcript,
	# File = "study.name_myPcalls.data.gedm.panp.txt"
		# P/A/M calls for all probe sets in each individual sample from the Presence/Absence Calling algorithm (PANP)
	# File = "study.name_myPcalls.data.gedm.panp.converted.to.binary.txt"
		# P/A/M calls for all probe sets in each individual sample from the Presence/Absence Calling algorithm (PANP) converted to binary calls i.e. P=1/A=0/M=0.5

	###if your arrayplatform is hgu133a or hgu133plus2 you will get an additional output file and plot (described below) for running the barcode presence/absence algorithm as a comparison to the PANP presence/absence algorithm implemented for processing the data
	# hgu133a Additional File = "study.name_myPcalls.data.gedm.barcode.GPL96.txt" 
	#	# P/A calls for all probe sets in each individual sample from the Presence/Absence Calling algorithm (barcode) for hgu133a microarray platform
	# hgu133plus2 Additional File = "study.name_myPcalls.data.gedm.barcode.GPL570.txt"
	#	# P/A calls for all probe sets in each individual sample from the Presence/Absence Calling algorithm (barcode) for hgu133plus2 microarray platform
	# Additional Plot = "study.name_PAcalls.distribution.PANP.vs.Barcode.png"
	#	# density plot showing distribution of Presence(P)/Absent(A)/Marginal(M) calls comparing the PANP algorithm to barcode algorithm produced only for the hgu133a and hgu133p microarray platforms
	###
	
	# File = "study.name_data.gedm.cases.controls.txt"
		# normalized gene expression data matrix (gedm) (rows = probes/probe sets/genes/features, columns = samples) for samples included in the demographic file and after removing outliers if any
	# File = "study.name_myPcalls.data.gedm.matrix.cases.controls.txt"
		# P/A/M calls for all probe sets in each individual sample from the Presence/Absence Calling algorithm (PANP) converted to binary calls i.e. P=1/A=0/M=0.5 for samples included in the demographic file and after removing outliers if any
	# File = "study.name_data.gedm.filtered.by.subjects.and.by.probesets.txt"
		# normalized gene expression data matrix (gedm) (rows = probes/probe sets/genes/features, columns = samples) filtered by subjects (outliers if any) and by probesets
	# File = "study.name_surrogate.variables.txt"
		# surrogate variables matrix for each individual
	# File = "study.name_residuals.data.gedm.post.sva.txt"
		# residual intensities matrix after removing effects of surrogate variables from the actual intensity i.e. (observed intensity - expected intensity) (rows = samples , columns = probes/probe sets/genes/features)
	# File = "study.name_residuals.data.gedm.transposed.post.sva.txt"
		# residual intensities matrix after removing effects of surrogate variables from the actual intensity i.e. (observed intensity - expected intensity) (rows = probes/probe sets/genes/features , columns = samples) - transposed
	# File = "study.name_probe.statistics.regression.residuals.post.sva.txt"
		#linear regression of residual intensities with groups for e.g. "1" (Group1 microarrays for e.g. "cases") and "0" (Group2 microarrays for e.g. "controls")
		# Columns-
		# standard error = standard_error_adjusted_res
		# beta estimate = coefficient_adjusted_res
		# fold change = fold_change_up_or_down_regulated_res
		# p-value = p_adjusted_res
	# File = "study.name_probe.statistics.regression.pre.and.post.sva.txt"
		#linear regression of actual intensities with groups for e.g. "1" (Group1 microarrays for e.g. "cases") and "0" (Group2 microarrays for e.g. "controls") with surrogate variables included as covariates in the model
		# Columns-
		# standard error unadjusted = standard_error_unadjusted
		# beta estimate unadjusted = coefficient_unadjusted
		# fold change unadjusted = fold_change_up_or_down_regulated_unadjusted
		# p-value unadjusted = p_unadjusted

		# standard error adjusted = standard_error_adjusted,
		# beta estimate adjusted = coefficient_adjusted,
		# fold change adjusted = fold_change_up_or_down_regulated_adjusted,
		# p-value adjusted = p_adjusted
	# File = "study.name_probe.statistics.regression.significants.pre.and.post.sva.p.less.than.0.05.txt"
		# output file above subset for probesets with association p-value adjusted < 0.05
	# Plot = "study.name_unadjusted.png"
		# 1) volcano plot 2) MA plot and 3) Histogram of association results pre-and post SVA
	# Plot = "study.name_adjusted.png"
		# 1) volcano plot 2) MA plot and 3) Histogram of association results pre-and post SVA
	# File = "study.name_residuals.data.gedm.transposed.post.sva.probesets.mapped.to.gene.txt"
		# residual intensities matrix after removing effects of surrogate variables from the actual intensity i.e. (observed intensity - expected intensity) (rows = probes/probe sets/genes/features , columns = samples) - transposed and mapped to gene
}

analyze = function(CELfile.directory.path = NULL, normalized.gedm.file.path = NULL, DemographicFileName.path = NULL, arrayplatform = NULL, outlier.samples = NULL, study.name = NULL, NSMP.no.transcript.file.path = NULL, NSMP.sense.transcript.file.path = NULL, NSMP.antisense.transcript.file.path = NULL, NSMP.overlap.transcript.file.path = NULL, PANP.script.path = NULL, probeset.to.gene.mapping.file.path = NULL){

if(CELfile.directory.path=="" && normalized.gedm.file.path==""){
	stop("\nAborting: PLease provide location of CEL files or name and location of gene expression data matrix\n\n")
}#end if

if(length(CELfile.directory.path)>0 && length(normalized.gedm.file)>0){
	stop("\nAborting: PLease provide location of CEL files or name and location of gene expression data matrix but not both\n\n")
}#end if

if(DemographicFileName.path==""){
	stop("\nAborting: PLease provide name and location of Demographic File\n\n")
}#end if

#if ((arrayplatform != "hgu133a") || (arrayplatform !="hgu133plus2") || (arrayplatform != "hgu95av2")) {
        #stop("\nAborting: arrayplatform type must be either hgu133a or hgu133plus2 or hgu95av2\n\n")
#}#end if

if(NSMP.no.transcript.file.path=="" || NSMP.sense.transcript.file.path=="" || NSMP.antisense.transcript.file.path=="" || NSMP.overlap.transcript.file.path=="")
{
	stop("\nAborting: PLease provide name and location of files with Negative Strand Matching Probesets (NSMP) identifiers identified as those matching to: NO TRANSCRIPTS, SENSE TRANSCRIPTS, ANTISENSE TRANSCRIPTS and OVERLAP TRANSCRIPTS\n\n")
	stop("\nThese files can be download from here: http://psychiatry.igm.jhmi.edu/wiki/index.php/Main_Page#Gene_Expression for hgu133a or hgu133plus2 or hgu95av2\n\n")
}#end if

if(PANP.script.path=="")
{
	stop("\nAborting: PLease provide name and location of script for Presence/Absence Calling (PANP)\n\n")
	stop("\nThis script can be download from here: http://psychiatry.igm.jhmi.edu/wiki/index.php/Main_Page#Gene_Expression\n\n")
}#end if


if(study.name==""){
	study.name = "study"
}#end if

	if(normalized.gedm.file.path==""){
		cat("\nReading: CEL files, processing into an Affy object, normalizing using Frozen Robust Multichip Analysis (fRMA) method\n\n")
		data = ReadAffy(celfile.path=CELfile.directory.path) #read in CEL files and create Affy object
		cat("\nSaving: Affy object for future use\n\n")
		save(data, file=paste(study.name,"_affy.RObject",sep=""))
		data #print details of data for e.g. microarray-platform, number of arrays/samples
		par(mfrow = c(1, 2),pty = "s")
		png(paste(study.name,"_boxplots.pre.and.post.normalization.png",sep=""))
		boxplot(data,cex=0,col="red",main=paste(study.name,"_pre.normalization",sep=""),las=2,cex.axis=0.70)
		data_frma = frma(data) #fRMA process data
		boxplot(data_frma,cex=0,col="red",main=paste(study.name,"_post.normalization",sep=""),las=2,cex.axis=0.70)
		dev.off()
		data_gedm = exprs(data_frma) #create gene expression data matrix (gedm), rows = probe sets, columns  arrays/samples
		write.table(data_gedm, file=paste(study.name,"_fRMA.normalized.txt",sep=""), sep="\t") #output gedm
		
	}#end if
	else{
		cat("\nReading: normalized gene expression data matrix\n\n")
		data_gedm = as.matrix(read.table(file=normalized.gedm.file.path, header=T, row.names=1, check.names=FALSE)) #read normalized gedm
		png(paste(study.name,"_boxplots.post.normalization.png",sep=""))
		boxplot(data_gedm,cex=0,col="red",main=paste(study.name,"_post.normalization",sep=""),las=2,cex.axis=0.70)
		dev.off()
	}#end else

		#present/absent/marginal probesets calls using PANP algorithm
		#references:
		#Warren P, Taylor D, Martini PGV, Jackson J, Bienkowska JR (Eds): Proceedings of the PANP - a New Method of Gene Detection on Oligonucleotide Expression Arrays
		#Proceedings of the 7th IEEE International Conference on Bioinformatics and Bioengineering, BIBE: October 14-17; Harvard Medical School, Boston, MA, USA. BIBE; 2007.
		cat("\nPresence/Absence Calling:",arrayplatform, "\n\n")
		png(filename=paste(study.name, "_NegativeStrandMatchingProbesets.intensity.distributions.png", sep=""))
	
		#plot the intensity distribution of 5 categories of Negative Strand Matching Probesets (NMSPs)
		plot(density(data_gedm), type='n',ylim=c(0,1), xlab="Log2(Intensity)", ylab="Probability density", main=paste(study.name, arrayplatform, "Distribution of NSMPs", sep="_")) 
		lines(density(data_gedm)) #ALL probesets distributions
	
		#extract list of Negative Strand Matching Probesets (NMSPs) identified by Warren P, Taylor D, Martini PGV, Jackson J, Bienkowska JR (Eds): Proceedings of the PANP - a New Method of Gene Detection on Oligonucleotide Expression Arrays
		#Proceedings of the 7th IEEE International Conference on Bioinformatics and Bioengineering, BIBE: October 14-17; Harvard Medical School, Boston, MA, USA. BIBE; 2007.
		#provided as part of the panp package in R
	
		if(arrayplatform == "hgu133a"){
			data(NSMPnames.hgu133a)
			NSMPnames.NSMP.PANP = data_gedm[NSMPnames.hgu133a,] #extract intensities for NSMPs identified by panp algorithm
			lines(density(NSMPnames.NSMP.PANP),col="red") #plot NSMPs identified by panp algorithm
			rm(NSMPnames.hgu133a, envir = globalenv()) #remove list of NSMPs identified by panp from the environment
		}#end if
		else if(arrayplatform == "hgu133plus2"){
			data(NSMPnames.hgu133plus2)
			NSMPnames.NSMP.PANP = data_gedm[NSMPnames.hgu133plus2,]
			lines(density(NSMPnames.NSMP.PANP),col="red") #plot NSMPs identified by panp algorithm
			rm(NSMPnames.hgu133plus2, envir = globalenv()) #remove list of NSMPs identified by panp from the environment
		}#end else if
		else {
		}#end else
		
		#Using Affymetrix annotation tables, we identified all probesets labeled as NSMPs that were not characterized with “cross hyb,” indicating the probeset may match to another gene.
		#We then mapped the NSMPs to the ENSEMBL transcript database version GRCh37.p8 and classified the NSMPs into four categories;
		#(i) probesets that did not map to a transcript at all;
		#(ii) probesets that detected sense transcripts;
		#(iii) probesets that detected antisense transcripts and
		#(iv) probesets that detected a sense transcript that overlaps with an antisense transcript
	
		#read list of 4 categories of Negative Strand Matching Probesets (NMSPs) provided as ADDITIONAL pre-processed files
		#define file path here
		NSMP.no.transcript = NSMP.no.transcript.file.path
		NSMP.sense.transcript =  NSMP.sense.transcript.file.path
		NSMP.antisense.transcript = NSMP.antisense.transcript.file.path
		NSMP.overlap.transcript = NSMP.overlap.transcript.file.path
	
		#(i) NSMPs that did not map to a transcript at all
		NSMPnames.NSMP.ensemble.mapped.no.transcript.table = as.matrix(read.table(file=NSMP.no.transcript,row.names=1))
		NSMPnames.NSMP.ensemble.mapped.no.transcript.names = dimnames(NSMPnames.NSMP.ensemble.mapped.no.transcript.table)[[1]]
		data.gedm.NSMP.ensemble.mapped.no.transcript = data_gedm[NSMPnames.NSMP.ensemble.mapped.no.transcript.names,]
		lines(density(data.gedm.NSMP.ensemble.mapped.no.transcript),col="green")
		
		#(ii) NSMPs that detected sense transcripts
		NSMPnames.NSMP.ensemble.mapped.sense.transcript.table = as.matrix(read.table(file=NSMP.sense.transcript,row.names=1))
		NSMPnames.NSMP.ensemble.mapped.sense.transcript.names = dimnames(NSMPnames.NSMP.ensemble.mapped.sense.transcript.table)[[1]]
		data.gedm.NSMP.ensemble.mapped.sense.transcript = data_gedm[NSMPnames.NSMP.ensemble.mapped.sense.transcript.names,]
		lines(density(data.gedm.NSMP.ensemble.mapped.sense.transcript),col="pink")
		
		#(iii) NSMPs that detected antisense transcripts
		NSMPnames.NSMP.ensemble.mapped.antisense.transcript.table = as.matrix(read.table(file=NSMP.antisense.transcript,row.names=1))
		NSMPnames.NSMP.ensemble.mapped.antisense.transcript.names = dimnames(NSMPnames.NSMP.ensemble.mapped.antisense.transcript.table)[[1]]
		data.gedm.NSMP.ensemble.mapped.antisense.transcript = data_gedm[NSMPnames.NSMP.ensemble.mapped.antisense.transcript.names,]
		lines(density(data.gedm.NSMP.ensemble.mapped.antisense.transcript),col="blue")
		
		#(iv) NSMPs that detected a sense transcript that overlaps with an antisense transcript
		NSMPnames.NSMP.ensemble.mapped.overlap.transcript.table = as.matrix(read.table(file=NSMP.overlap.transcript,row.names=1))
		NSMPnames.NSMP.ensemble.mapped.overlap.transcript.names = dimnames(NSMPnames.NSMP.ensemble.mapped.overlap.transcript.table)[[1]]
		data.gedm.NSMP.ensemble.mapped.overlap.transcript = data_gedm[NSMPnames.NSMP.ensemble.mapped.overlap.transcript.names,]
		lines(density(data.gedm.NSMP.ensemble.mapped.overlap.transcript),col="orange")
		
		legend("topright", legend=c("all_probesets","PANP_NSMPs","ENSEMBLE_NSMPs_no_transcript","ENSEMBLE_NSMPs_sense_transcript","ENSEMBLE_NSMPs_antisense_transcript","ENSEMBLE_NSMPs_overlap_transcript"),cex=0.8,col=c('black','red','green','pink','blue','orange'),lty=1,title="NSMPs")
		
		dev.off()
	
		#run panp algorithm using #(i) NSMPs that did not map to a transcript at all
		source(PANP.script.path)
		NSMPnames = NSMPnames.NSMP.ensemble.mapped.no.transcript.names #use NSMPs list defined by our algorithm
		data_gedm_panp = pa.calls.mod(data_gedm,NSMP.names=NSMPnames, chip=arrayplatform) #make present/absent calls using panp
		myPcalls_data_gedm_panp = data_gedm_panp$Pcalls #extract present/absent calls
		write.table(myPcalls_data_gedm_panp,file=paste(study.name,"_myPcalls.data.gedm.panp.txt",sep="")) #write present/absent calls to a file
	
		if(arrayplatform == "hgu133a"){
			#compare panp present/absent calls against barcode algorithm present/absent calls
			data_gedm_barcode = barcode(data_gedm, output="binary", platform="GPL96")
			write.table(data_gedm_barcode,file=paste(study.name,"_myPcalls.data.gedm.barcode.GPL96.txt",sep=""))
			
			#convert present/absent calls from panp to binary
			myPcalls_data_gedm_panp_reload = as.matrix(read.table(file=paste(study.name,"_myPcalls.data.gedm.panp.txt",sep=""),header=T,row.names=1, check.names=FALSE))
			myPcalls_data_gedm_panp_reload[myPcalls_data_gedm_panp_reload=="A"]=0 #absent
			myPcalls_data_gedm_panp_reload[myPcalls_data_gedm_panp_reload=="P"]=1 #present
			myPcalls_data_gedm_panp_reload[myPcalls_data_gedm_panp_reload=="M"]=0.5 #marginal
			myPcalls_data_gedm_panp_reload_as_numeric = apply(myPcalls_data_gedm_panp_reload,c(1,2),as.numeric) #convert to numeric matrix
			write.table(myPcalls_data_gedm_panp_reload_as_numeric,file=paste(study.name,"_myPcalls.data.gedm.panp.converted.to.binary.txt",sep="")) #write numeric present/absent calls to a file
	
			#plot density of present absent calls from panp vs. barcode comparison
			png(filename=paste(study.name,"_PAcalls.distribution.PANP.vs.Barcode.png",sep=""))
		
			plot(density(data_gedm_barcode), type='n', xlab="present=1,absent=0,marginal=0.5", main=paste(study.name,arrayplatform,"Dist. of PAcalls_PANP_vs_Barcode",sep="_"), cex.main=0.90)
			lines(density(data_gedm_barcode))
			lines(density(myPcalls_data_gedm_panp_reload_as_numeric),col='red')
				
			legend("topright", legend=c("Barcode_PAcalls","PANP_PAcalls"),cex=0.8,col=c('black','red'),lty=1,title="PAcalls")
		
			dev.off()
		}#end if
		else if(arrayplatform == "hgu133plus2"){
			#compare panp present/absent calls against barcode algorithm present/absent calls
			data_gedm_barcode = barcode(data_gedm, output="binary", platform="GPL570")
			write.table(data_gedm_barcode,file=paste(study.name,"_myPcalls.data.gedm.barcode.GPL570.txt",sep=""))

			#convert present/absent calls from panp to binary
			myPcalls_data_gedm_panp_reload = as.matrix(read.table(file=paste(study.name,"_myPcalls.data.gedm.panp.txt",sep=""),header=T,row.names=1))
			myPcalls_data_gedm_panp_reload[myPcalls_data_gedm_panp_reload=="A"]=0 #absent
			myPcalls_data_gedm_panp_reload[myPcalls_data_gedm_panp_reload=="P"]=1 #present
			myPcalls_data_gedm_panp_reload[myPcalls_data_gedm_panp_reload=="M"]=0.5 #marginal
			myPcalls_data_gedm_panp_reload_as_numeric = apply(myPcalls_data_gedm_panp_reload,c(1,2),as.numeric) #convert to numeric matrix
			write.table(myPcalls_data_gedm_panp_reload_as_numeric,file=paste(study.name,"_myPcalls.data.gedm.panp.converted.to.binary.txt",sep="")) #write numeric present/absent calls to a file
	
			#plot density of present absent calls from panp vs. barcode comparison
			png(filename=paste(study.name,"_PAcalls.distribution.PANP.vs.Barcode.png",sep=""))
		
			plot(density(data_gedm_barcode), type='n', xlab="present=1,absent=0,marginal=0.5", main=paste(study.name,arrayplatform,"Dist. of PAcalls_PANP_vs_Barcode",sep="_"), cex.main=0.90)
			lines(density(data_gedm_barcode))
			lines(density(myPcalls_data_gedm_panp_reload_as_numeric),col='red')
				
			legend("topright", legend=c("Barcode_PAcalls","PANP_PAcalls"),cex=0.8,col=c('black','red'),lty=1,title="PAcalls")
		
			dev.off()

		}#end else if
		else {
		}#end else
	
		cat("\nAnalyzing: normalized gene expression data matrix\n\n")

		#START ANALYSIS (sva, differential expression)

		#An m x 2 matrix with m individuals/arrays divided into groups. 
		#Column1 = Unique Identifier (ID) matching the column header in your gene expression data matrix (gedm)
		#Column2 = 1 (Group1 arrays for e.g. "cases") and 0 (Group2 arrays for e.g. "controls")
		Demographic.matrix = as.matrix(read.table(file=DemographicFileName.path,row.names=1,check.names=FALSE))

		data_gedm_sample_IDs = dimnames(data_gedm)[[2]]	#sample_IDs (header row in your gene expression data matrix)

		#remove outliers (if any, IDs should match sample_IDs (header row in your gene expression data matrix))
		if(length(outlier.samples)>0){
			cat("\nAnalyzing: removing outliers from Demographic File\n\n")
			check_outlier_ids = all(outlier.samples %in% data_gedm_sample_IDs)
			if(check_outlier_ids==TRUE){
				t = as.matrix(Demographic.matrix[!rownames(Demographic.matrix) %in% outlier.samples,])
				Demographic.matrix = t
			}#end if	
			else{
				stop("\nAborting: Your subject IDs do not match between your outlier ID list and gene expression data matrix. No outliers excluded\n\n")
			}#end else
			
		}#end if		

		Demographic.matrix.cases.controls = dimnames(Demographic.matrix)[[1]] #sample_IDs (first column in the Demographic File)

		#check to see if samples in Demographic file exist in the gedm
		check_ids = all(Demographic.matrix.cases.controls %in% data_gedm_sample_IDs)
		if(check_ids==TRUE){
			cat("\nAnalyzing: subsetting gene expression data matrix using Demographic File\n\n")
			data_gedm.cases.controls = subset(data_gedm,select=Demographic.matrix.cases.controls)
			s = data_gedm.cases.controls
			write.table(data_gedm.cases.controls,file=paste(study.name,"_data.gedm.cases.controls.txt",sep=""))
		}#end if
		else{
			stop("\nAborting: Your subject IDs do not match between your demographic file and gene expression data matrix. Please check your input files\n\n")
		}#end else

		cat("\nAnalyzing: filter probesets using Presence/Absence calls\n\n")
		#read in present/marginal/absent calls file generated using PANP algorithm (arrayplatform specific above) for all probesets and individuals
		data.gedm.PANP.calls = as.matrix(read.table(file=paste(study.name,"_myPcalls.data.gedm.panp.converted.to.binary.txt",sep=""),header=T,row.names=1, check.names=FALSE))
		data.gedm.PANP.calls.cases.controls = subset(data.gedm.PANP.calls,select=Demographic.matrix.cases.controls)
		g = data.gedm.PANP.calls.cases.controls
		write.table(data.gedm.PANP.calls.cases.controls,file=paste(study.name,"_myPcalls.data.gedm.matrix.cases.controls.txt",sep=""))
	
		#create group assignment variable (cases = 1 and controls = 0)
		grouplabels = as.vector(Demographic.matrix[,1])
		grp = as.numeric(grouplabels)

		cases = length(grp[grp==1])
	
		#filter by present, absent & marginal calls
		present_cases = apply(g, 1, function(x) (sum(x[1:cases]))) #calculate sum of present calls in cases for each probeset
		present_controls = apply(g, 1, function(x) (sum(x[(cases+1):length(x)]))) #calculate sum of present calls in controls for each probeset
	
		j=cbind(s,present_cases,present_controls) #combine sum of calls with gedm
		k=subset(j, present_cases>0 | present_controls>0) #keep probesets with a single present call across all samples
		l=subset(k, select = -c(present_cases,present_controls)) #removing the last 2 columns (sum of PA calls in cases and sum of PA calls in controls)
		s=l
	
		write.table(s,file=paste(study.name,"_data.gedm.filtered.by.subjects.and.by.probesets.txt",sep=""), sep = "\t") 
	
		cat("\nAnalyzing: association (linear regression) pre-surrogate variable analysis\n\n")	

		#linear regression (association BP_CONTROL)
		stats = apply(s, 1, function(x) linear_regression(x,grp))
		tstats = t(stats)
		p_unadjusted = tstats[,2]
		coefficient_unadjusted = tstats[,1]
		standard_error_unadjusted = tstats[,3]
		
		fold_change_unadjusted = as.vector(2^coefficient_unadjusted) #calculate actual FC
		h = cbind(fold_change_unadjusted)
		
		fold_change_up_or_down_regulated_unadjusted = apply(h, 1, function(x) if(x[1]<1) x[1]=-(1/x[1]) else x[1]=x[1]) #change sign of FC, from above: if(coefficient is negative, FC<1) if(coefficient is positive, FC>1), now, if(FC<1) take inverse and put a negative sign infront of FC, esle leave FC as is
		
		cat("\nAnalyzing: surrogate variable analysis (SVA)\n\n")	

		#surrogate variable analysis
		mod = cbind(rep(1,length(grp)),grp)
		mod0 = cbind(rep(1,length(grp)))
		svaobj = sva(s,mod,mod0,method="irw")
	
		write.table(svaobj$sv,file=paste(study.name,"_surrogate.variables.txt",sep=""), sep = "\t") 
	
		cat("\n\nAnalyzing: linear regression of surrogate variables to build residual matrix\n\n")

		#linear_regression_adjusted (BP_CONTROL) to build residual matrix
		s_trans = t(s)
		res = apply(s_trans, 2, function(x) linear_regression_adjusted_for_residuals(x,svaobj$sv))
		res_matrix = cbind(res)
		res_matrix_w_grp_studyid = cbind(grp,study.name,res_matrix)
		res_matrix_w_grp_study_id_df = data.frame(res_matrix_w_grp_studyid, check.names=FALSE)
	
		write.table(res_matrix_w_grp_study_id_df,file=paste(study.name,"_residuals.data.gedm.post.sva.txt",sep=""), sep = "\t")
		write.table(t(res_matrix_w_grp_study_id_df),file=paste(study.name,"_residuals.data.gedm.transposed.post.sva.txt",sep=""), sep = "\t")
	
		cat("\nAnalyzing: association (linear regression) post-surrogate variable analysis with residual matrix\n\n")	
	
		#linear regression (association BP_CONTROL with residual matrix)
		stats_res = apply(res_matrix, 2, function(x) linear_regression(x,grp))
		tstats_res = t(stats_res)
		p_adjusted_res = tstats_res[,2]
		coefficient_adjusted_res = tstats_res[,1]
		standard_error_adjusted_res = tstats_res[,3]
	
		fold_change_adjusted_res = as.vector(2^coefficient_adjusted_res) #calculate actual FC
		h_res = cbind(fold_change_adjusted_res)
	
		fold_change_up_or_down_regulated_res = apply(h_res, 1, function(x) if(x[1]<1) x[1]=-(1/x[1]) else x[1]=x[1]) #change sign of FC, from above: if(coefficient is negative, FC<1) if(coefficient is positive, FC>1), now, if(FC<1) take inverse and put a negative sign infront of FC, esle leave FC as is
	
		foo = cbind(standard_error_adjusted_res,coefficient_adjusted_res,fold_change_up_or_down_regulated_res,p_adjusted_res)
		write.table(foo,file=paste(study.name,"_probe.statistics.regression.residuals.post.sva.txt",sep=""), sep="\t")
	
		cat("\nAnalyzing: association (linear regression) post-surrogate variable adjusting for surrogate variables as covariates\n\n")	

		#linear_regression_adjusted (assaciation BP_CONTROL)
		stats_a = apply(s, 1, function(x) linear_regression_adjusted(x,grp,svaobj$sv))
		tstats_a = t(stats_a)
		p_adjusted = tstats_a[,2]
		coefficient_adjusted = tstats_a[,1]
		standard_error_adjusted = tstats_a[,3]
	
		fold_change_adjusted = as.vector(2^coefficient_adjusted) #calculate actual FC
		h = cbind(fold_change_adjusted)
		
		fold_change_up_or_down_regulated_adjusted = apply(h, 1, function(x) if(x[1]<1) x[1]=-(1/x[1]) else x[1]=x[1]) #change sign of FC, from above: if(coefficient is negative, FC<1) if(coefficient is positive, FC>1), now, if(FC<1) take inverse and put a negative sign infront of FC, esle leave FC as is
	
		e = cbind(standard_error_unadjusted,coefficient_unadjusted,fold_change_up_or_down_regulated_unadjusted,p_unadjusted,standard_error_adjusted,coefficient_adjusted,fold_change_up_or_down_regulated_adjusted,p_adjusted)
		edf = data.frame(e)
	
		e_significants = subset(edf, p_adjusted<0.05)
	
		write.table(e,file=paste(study.name,"_probe.statistics.regression.pre.and.post.sva.txt",sep=""), sep = "\t")
	
		write.table(e_significants,file=paste(study.name,"_probe.statistics.regression.significants.pre.and.post.sva.p.less.than.0.05.txt",sep=""), sep = "\t")
	
		cat("\nAnalyzing: generating 1) volcano plot 2) MA plot and 3) Histogram of association results pre-and post SVA \n\n")	

		p_unadjusted.trans = -1 * log10(p_unadjusted) 
		p_adjusted.trans = -1 * log10(p_adjusted) 
	
		M=coefficient_unadjusted
		M_a=coefficient_adjusted
	
		png(filename=paste(study.name,"_unadjusted.png",sep=""), width = 1400, height = 800)
		par(mfrow = c(1, 3),pty = "s")
	
		#volcano plot
		plot(range(M),range(p_unadjusted.trans),type="n",ylab="-log10(p-value_unadjusted)",xlab="Beta Estimate",main=paste(study.name,"_Unadjusted Volcano Plot",sep="")) 
		points(M,p_unadjusted.trans,col="black",cex=0.2, pch=16)
		points(M[(p_unadjusted.trans > 1.3 & M > 1)],p_unadjusted.trans[(p_unadjusted.trans > 1.3 & M > 1)],col="red",pch=16) 
		points(M[(p_unadjusted.trans > 1.3 & M < -1)],p_unadjusted.trans[(p_unadjusted.trans > 1.3 & M < -1)],col="green",pch=16) 
		abline(h = 1.3)
		abline(v = -1)
		abline(v = 1)
	
		#ma plot variables
		cases = length(grp[grp==1]) 
		controls = length(grp[grp==0]) 
		cases.m = apply(s[,1:cases],1,mean)
		controls.m = apply(s[,controls:length(grp)],1,mean)
	
		#ma plot
		a = (cases.m + controls.m)/2.0
		plot(range(a),range(coefficient_unadjusted),type="n",ylab="M",xlab="A",main=paste(study.name,"_Unadjusted MA Plot",sep=""))
		points(a,coefficient_unadjusted,col="black",cex=0.2, pch=16)
		points(a[(coefficient_unadjusted > 1)],coefficient_unadjusted[(coefficient_unadjusted > 1)],col="red",pch=16) 
		points(a[(coefficient_unadjusted < -1)],coefficient_unadjusted[(coefficient_unadjusted < -1)],col="green",pch=16) 
		abline(h = -1)
		abline(h = 1)
	
		hist(p_unadjusted, main = paste(study.name,"_Unajusted P-values",sep=""), xlab = "P-value", col = "grey")
		
		dev.off()
	
		png(filename=paste(study.name,"_adjusted.png",sep=""), width = 1400, height = 800)
	
		par(mfrow = c(1, 3),pty = "s")
	
		#volcano plot
		plot(range(M_a),range(p_unadjusted.trans),type="n",ylab="-log10(p-value_adjusted)",xlab="M",main=paste(study.name,"_adjusted Volcano Plot",sep="")) 
		points(M_a,p_adjusted.trans,col="black",cex=0.2, pch=16)
		points(M_a[(p_adjusted.trans > 1.3 & M_a > 1)],p_adjusted.trans[(p_adjusted.trans > 1.3 & M_a > 1)],col="red",pch=16) 
		points(M_a[(p_adjusted.trans > 1.3 & M_a < -1)],p_adjusted.trans[(p_adjusted.trans > 1.3 & M_a < -1)],col="green",pch=16) 
		abline(h = 1.3)
		abline(v = -1)
		abline(v = 1)
	
		#ma plot
		a = (cases.m + controls.m)/2.0
		plot(range(a),range(coefficient_adjusted),type="n",ylab="M",xlab="A",main=paste(study.name,"_Adjusted MA Plot",sep="")) 
		points(a,coefficient_adjusted,col="black",cex=0.2, pch=16)
		points(a[(coefficient_adjusted > 1)],coefficient_adjusted[(coefficient_adjusted > 1)],col="red",pch=16) 
		points(a[(coefficient_adjusted < -1)],coefficient_adjusted[(coefficient_adjusted < -1)],col="green",pch=16) 
		abline(h = -1)
		abline(h = 1)
	
		hist(p_adjusted, main = paste(study.name,"_Adjusted P-values",sep=""), xlab = "P-value", col = "grey")
	
		dev.off()

		cat("\nAnalyzing: mapping probesets to gene for residual matrix which can be used for association or mega-analysis\n\n")
		probesets_to_gene_mapping_file = as.matrix(read.table(file=probeset.to.gene.mapping.file.path, header=T, row.names=1, check.names=FALSE)) #read probeset to gene mapping
		data_gedm_residuals_merged_with_probeset_to_gene_mapping_file = merge(probesets_to_gene_mapping_file,t(res_matrix_w_grp_study_id_df),by=c(0,0))
		data_gedm_residuals_merged_with_probeset_to_gene_mapping_file_remove_some_columns=subset(data_gedm_residuals_merged_with_probeset_to_gene_mapping_file, select = -c(Row.names))
		write.table(data_gedm_residuals_merged_with_probeset_to_gene_mapping_file_remove_some_columns,file=paste(study.name,"_residuals.data.gedm.transposed.post.sva.probesets.mapped.to.gene.txt",sep=""), sep="\t")
		
}#end function analyze

# linear regression analysis for data processing
linear_regression = function(x,grp){
	lm_model = lm(x~grp)
	coefficient = coefficients(summary(lm_model))[,1][2]
	unadjusted_P = coefficients(summary(lm_model))[,4][2]
	se = coefficients(summary(lm_model))[,2][2]
	unadjusted_stats = cbind(coefficient,unadjusted_P,se)
}

linear_regression_adjusted = function(x,grp,sv){
	adjusted_lm_model = lm(x~grp+sv)
	adjusted_coefficient = coefficients(summary(adjusted_lm_model))[,1][2]
	adjusted_P = coefficients(summary(adjusted_lm_model))[,4][2]
	adjusted_se = se = coefficients(summary(adjusted_lm_model))[,2][2]
	adjusted_stats = cbind(adjusted_coefficient,adjusted_P,adjusted_se)
}

linear_regression_adjusted_for_residuals = function(x,sv){
	adjusted_lm_model_for_residuals = lm(x~sv)
	residuals = adjusted_lm_model_for_residuals$residuals
}

#execute functions
load_required_packages()

variables = define_paths_and_variables()

analyze(CELfile.directory.path = NULL, normalized.gedm.file.path = variables$normalized.gedm.file.path, DemographicFileName.path = variables$DemographicFileName.path, arrayplatform = variables$arrayplatform, outlier.samples = as.vector(variables$outlier.samples), study.name = variables$study.name, NSMP.no.transcript.file.path = variables$NSMP.no.transcript.file.path, NSMP.sense.transcript.file.path = variables$NSMP.sense.transcript.file.path, NSMP.antisense.transcript.file.path = variables$NSMP.antisense.transcript.file.path, NSMP.overlap.transcript.file.path = variables$NSMP.overlap.transcript.file.path, PANP.script.path = variables$PANP.script.path, probeset.to.gene.mapping.file.path = variables$probeset.to.gene.mapping.file.path)







