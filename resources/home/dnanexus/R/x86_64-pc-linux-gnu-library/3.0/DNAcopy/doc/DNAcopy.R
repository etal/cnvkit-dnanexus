### R code from vignette source 'DNAcopy.Rnw'

###################################################
### code chunk number 1: DNAcopy.Rnw:74-75
###################################################
library(DNAcopy)


###################################################
### code chunk number 2: DNAcopy.Rnw:78-79
###################################################
data(coriell)


###################################################
### code chunk number 3: DNAcopy.Rnw:85-88
###################################################
CNA.object <- CNA(cbind(coriell$Coriell.05296),
                  coriell$Chromosome,coriell$Position,
                  data.type="logratio",sampleid="c05296")


###################################################
### code chunk number 4: DNAcopy.Rnw:96-97
###################################################
smoothed.CNA.object <- smooth.CNA(CNA.object)


###################################################
### code chunk number 5: DNAcopy.Rnw:105-106
###################################################
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)


###################################################
### code chunk number 6: DNAcopy.Rnw:120-121
###################################################
plot(segment.smoothed.CNA.object, plot.type="w")


###################################################
### code chunk number 7: DNAcopy.Rnw:129-130
###################################################
plot(segment.smoothed.CNA.object, plot.type="s") 


###################################################
### code chunk number 8: DNAcopy.Rnw:157-158
###################################################
plot(segment.smoothed.CNA.object, plot.type="p")


###################################################
### code chunk number 9: DNAcopy.Rnw:169-172
###################################################
sdundo.CNA.object <- segment(smoothed.CNA.object, 
                             undo.splits="sdundo", 
                             undo.SD=3,verbose=1)


###################################################
### code chunk number 10: DNAcopy.Rnw:177-178
###################################################
plot(sdundo.CNA.object,plot.type="s")


