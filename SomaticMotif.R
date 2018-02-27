.libPaths("/gpfs/users/yanghao/software/anaconda2/lib/R/library")
suppressPackageStartupMessages(library("deconstructSigs"))

Arg <- commandArgs();
nameforinput = Arg[6]
annotationforinput = Arg[7]
Command = paste("/gpfs/users/yanghao/project/somatic_motif_from_bam/SomaticMotif.sh", nameforinput, sep=" ")
system(Command)

fileforinput = paste(nameforinput, ".somatic.forMotif",sep="")
dat = read.delim(fileforinput,header=TRUE)

sigs.input <- mut.to.sigs.input(mut.ref = dat, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt")

SampleSig = whichSignatures(tumor.ref = sigs.input, 
                           signatures.ref = signatures.cosmic, 
                           sample.id = nameforinput, 
                           contexts.needed = TRUE,
                           tri.counts.method = 'exome')

outputname = paste(nameforinput, ".SomaticMotif.pdf",sep="")
pdf(outputname)
plotSignatures(SampleSig, sub = annotationforinput)
makePie(SampleSig)
dev.off()
