# install and load necessary libraries

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rtracklayer")

library(rtracklayer)

# 4.1. Set you working directory:
# set you working directory (you can also use the '...' under the 'File' tab on 
# the right of Rstudio to navigate to your folder then click 'More' > 'Set as 
# working directory')

setwd("G:/My Drive/Workshop_pan_genomics_12042023")

# 4.2. Load the Minimap2 output

miniMap_out <- read.table('Minimap2_out/Minimap_melp_erato.sam', header = FALSE, sep = '\t')

# set the column names
colnames(miniMap_out) <- c('queryName', 'queryLength', 'queryStart', 'queryEnd', 
                           'char', 
                           'targetName', 'targetLength', 'targetStart', 'targetEnd', 
                           'matchingBases', 'matchLength', 'matchQuality')

# check the table
head(miniMap_out)
  
# 4.3. Plot the first Minimap2 match

# First define x-axis coordinates. This will help us later define the genomic window (xlim) we want to zoom in on.
start = 0
end = 3000000

# Plot an empty plot (so we can fill it with rectangles (~genome) and polygons (~alignments)).
# For the y-axis (ylim) coordinates I arbitrarily use c(0,10).
plot(NULL, xlim = c(start, end), ylim = c(0,10), axes=F, ylab = '', xlab = '')

# Now we can draw two simple rectangles, one will define the genomic interval/sequence of H. melpomene (deepskyblue), the other H. erato (mediumseagreen).
rect(0,8,2000000,9, col = 'deepskyblue', border = NA)
rect(0,1,2000000,2, col = 'mediumseagreen', border = NA)

# We also know the positions of the optix gene and we can add these with the same rect() function trick.
rect(705604,9,706407,10, col = 'red', border = 'red') # position optix melpomene (only has one exon)

rect(1239943,0,1239972,1, col = 'red', border = 'red') # first exon position optix erato
rect(1250591,0,1251211,1, col = 'red', border = 'red') # second exon position optix erato
rect(1239943,0.5,1251211,0.5, col = 'red', border = 'red') # A little line between the two and we have a gene model!

# With the text function, we can add the gene and species names at the appropriate coordinates.
text(706407, 9.5, substitute(paste(italic('optix'))), pos = 4)
text(1251211, 9.5, substitute(paste(italic('optix'))), pos = 4)

text(start, 7.5, substitute(paste(italic('H. melpomene'))), pos = 4)
text(start, 2.5, substitute(paste(italic('H. erato'))), pos = 4)

# Finally, with the polygon function, we map the alignment connections between the sequences of H. melpomene (target) and H. erato (query).
polygon(x = c(miniMap_out$targetStart[1], miniMap_out$targetEnd[1], miniMap_out$queryEnd[1], miniMap_out$queryStart[1]), 
        y = c(8,8,2,2),
        col = adjustcolor('black', alpha.f = 0.1), border = FALSE)
###

# 4.4. Plot all Minimap2 match

# First define x-axis coordinates. This will help us later define the genomic window (xlim) we want to zoom in on.
start = 0
end = 3000000

# Plot an empty plot (so we can fill it with rectangles (~genome) and polygons (~alignments)).
# For the y-axis (ylim) coordinates I arbitrarily use c(0,10).
plot(NULL, xlim = c(start, end), ylim = c(0,10), axes=F, ylab = '', xlab = '')

# Now we can draw two simple rectangles, one will define the genomic interval/sequence of H. melpomene (deepskyblue), the other H. erato (mediumseagreen).
rect(0,8,2000000,9, col = 'deepskyblue', border = NA)
rect(0,1,2000000,2, col = 'mediumseagreen', border = NA)

# We also know the positions of the optix gene and we can add these with the same rect() function trick.
rect(705604,9,706407,10, col = 'red', border = 'red') # position optix melpomene (only has one exon)

rect(1239943,0,1239972,1, col = 'red', border = 'red') # first exon position optix erato
rect(1250591,0,1251211,1, col = 'red', border = 'red') # second exon position optix erato
rect(1239943,0.5,1251211,0.5, col = 'red', border = 'red') # A little line between the two and we have a gene model!

# With the text function, we can add the gene and species names at the appropriate coordinates.
text(706407, 9.5, substitute(paste(italic('optix'))), pos = 4)
text(1251211, 9.5, substitute(paste(italic('optix'))), pos = 4)

text(start, 7.5, substitute(paste(italic('H. melpomene'))), pos = 4)
text(start, 2.5, substitute(paste(italic('H. erato'))), pos = 4)

# Loop through each row in the Minimap2 alignment table and plot each match as a polygon.
for(e in 1:nrow(miniMap_out)){
  
    polygon(x = c(miniMap_out$targetStart[e], miniMap_out$targetEnd[e], miniMap_out$queryEnd[e], miniMap_out$queryStart[e]), 
            y = c(8,8,2,2),
            col = adjustcolor('black', alpha.f = 0.1), border = FALSE)
}
###


# 4.5. Modify relative alignment
# Let's make sure optix is aligned in the center of the plot between the two species

# First define x-axis coordinates. This will help us later define the genomic window (xlim) we want to zoom in on.
start = 0
end = 3000000

# Calculate the offset between the positions of optix in the two fasta sequences.
# We will use this to modify all the x-axes coordinates for H. melpomene.
plotDiff = 1251211 - 706407

# Plot an empty plot (so we can fill it with rectangles (~genome) and polygons (~alignments)).
# For the y-axis (ylim) coordinates I arbitrarily use c(0,10).
plot(NULL, xlim = c(start, end), ylim = c(0,10), axes=F, ylab = '', xlab = '')

# Now we can draw two simple rectangles, one will define the genomic interval/sequence of H. melpomene (deepskyblue), the other H. erato (mediumseagreen).
rect(0+plotDiff,8,2000000+plotDiff,9, col = 'deepskyblue', border = NA)
rect(0,1,2000000,2, col = 'mediumseagreen', border = NA)

# We also know the positions of the optix gene and we can add these with the same rect() function trick.
rect(705604+plotDiff,9,706407+plotDiff,10, col = 'red', border = 'red') # position optix melpomene (only has one exon)

rect(1239943,0,1239972,1, col = 'red', border = 'red') # first exon position optix erato
rect(1250591,0,1251211,1, col = 'red', border = 'red') # second exon position optix erato
rect(1239943,0.5,1251211,0.5, col = 'red', border = 'red') # A little line between the two and we have a gene model!

# With the text function, we can add the gene and species names at the appropriate coordinates.
text(706407+plotDiff, 9.5, substitute(paste(italic('optix'))), pos = 4)
text(1251211, 9.5, substitute(paste(italic('optix'))), pos = 4)

text(start, 7.5, substitute(paste(italic('H. melpomene'))), pos = 4)
text(start, 2.5, substitute(paste(italic('H. erato'))), pos = 4)

# Plot the polygons representing the alignments.
for(e in 1:nrow(miniMap_out)){
  
  polygon(x = c(miniMap_out$targetStart[e]+plotDiff, miniMap_out$targetEnd[e]+plotDiff, miniMap_out$queryEnd[e], miniMap_out$queryStart[e]), 
          y = c(8,8,2,2),
          col = adjustcolor('black', alpha.f = 0.1), border = FALSE)
}
###


# 4.6. Zoom in 
###
# Exercise: Zoom in on 10,000 bp before and 200,000 bp after optix

# but before we do this, let's build the layout of our plot first

layout(matrix(c(1,4,5,3,6,7,2), nrow=7, byrow=TRUE), height = c(0.5,1,1,2,1,1,0.5))
layout.show(n=7)

par(mar = c(0.5,5,0.5,1), xpd = FALSE)

# Define the interval (relative to the H. erato genome)
start = 1239943 - 10000
end = 1251211 + 200000

# Calculate the offset between the positions of optix in the two fasta sequences.
# We will use this to modify all the x-axes coordinates for H. melpomene.
plotDiff = 1251211 - 706407

# Plot an empty plot (so we can fill it with rectangles (~genome) and polygons (~alignments)).
# For the y-axis (ylim) coordinates I arbitrarily use c(0,10).
plot(NULL, xlim = c(start, end), ylim = c(0,10), axes=F, ylab = '', xlab = '')

# Now we can draw two simple rectangles, one will define the genomic interval/sequence of H. melpomene (deepskyblue), the other H. erato (mediumseagreen).
rect(0+plotDiff,8,2000000+plotDiff,9, col = 'deepskyblue', border = NA)
rect(0,1,2000000,2, col = 'mediumseagreen', border = NA)

# We also know the positions of the optix gene and we can add these with the same rect() function trick.
rect(705604+plotDiff,9,706407+plotDiff,10, col = 'red', border = 'red') # position optix melpomene (only has one exon)

rect(1239943,0,1239972,1, col = 'red', border = 'red') # first exon position optix erato
rect(1250591,0,1251211,1, col = 'red', border = 'red') # second exon position optix erato
rect(1239943,0.5,1251211,0.5, col = 'red', border = 'red') # A little line between the two and we have a gene model!

# With the text function, we can add the gene and species names at the appropriate coordinates.
text(706407+plotDiff, 9.5, substitute(paste(italic('optix'))), pos = 4)
text(1251211, 9.5, substitute(paste(italic('optix'))), pos = 4)

text(start, 7.5, substitute(paste(italic('H. melpomene'))), pos = 4)
text(start, 2.5, substitute(paste(italic('H. erato'))), pos = 4)

for(e in 1:nrow(miniMap_out)){
  
  # I add here a filter so that alignments outside of the plotting area don't get included.
  if(miniMap_out$targetStart[e]+plotDiff > start & miniMap_out$targetEnd[e]+plotDiff < end & miniMap_out$queryStart[e] > start & miniMap_out$queryEnd[e] < end){
  
    polygon(x = c(miniMap_out$targetStart[e]+plotDiff, miniMap_out$targetEnd[e]+plotDiff, miniMap_out$queryEnd[e], miniMap_out$queryStart[e]), 
          y = c(8,8,2,2),
          col = adjustcolor('black', alpha.f = miniMap_out$matchingBases[e]/miniMap_out$matchLength[e]), border = FALSE)
    
  }
}

mtext('Minimap2 alignment', side = 2, cex=0.8, padj = -1, col = 'black')
###



# 4.7. Add transposable element annotations
# We can plot Transposable Elements on top of this

# Load the TE tables (only the columns that are relevant to us).
TE_erato <- read.table('TEs/H_erato_1801_TE.txt', header = FALSE)[,c(2,5:7)]
TE_melp <- read.table('TEs/H_melp_18003_TE.txt', header = FALSE)[,c(2,5:7)]

# Fix the column names
colnames(TE_erato) <- c('subs','scaf', 'start', 'end')
colnames(TE_melp) <- c('subs','scaf', 'start', 'end')

# (optional) remove TEs with large number of substitutions compared to database
TE_erato <- subset(TE_erato, TE_erato$subs < 15)
TE_melp  <- subset(TE_melp, TE_melp$subs < 15)

# Loop through the rows and plot each TE as a small orange rectangle.

for(e in 1:nrow(TE_melp)){
  rect(TE_melp$start[e]+plotDiff,8.8,TE_melp$end[e]+plotDiff,9, col = 'orange', border = NA)
}

for(e in 1:nrow(TE_erato)){
  rect(TE_erato$start[e],1,TE_erato$end[e],1.2, col = 'orange', border = NA)
}

# Add a little text note on the plot
mtext('TE', side = 1, cex=0.8, padj = -3, las = 1, adj=1, col = 'orange')


# 4.8. Add ATAC-seq data tracks
# Add the ATAC-seq data of some butterfly brain and wing tissue

# Load the ATAC-seq data (make sure the 'rtracklayer' is loaded).

erato_5th_brain <- import.bedGraph("ATAC/brain_5th_H_erato_normalized_mean.w30s0bin.bg")
erato_5th_FW <- import.bedGraph("ATAC/FW_5th_H_erato_normalized_mean.w30s0bin.bg")

melp_5th_brain <- import.bedGraph("ATAC/brain_5th_H_melp_normalized_mean.w30s0bin.bg")
melp_5th_FW <- import.bedGraph("ATAC/FW_5th_H_melp_normalized_mean.w30s0bin.bg")

# Sample 1

# Plot an empty plot (remember, this will fill a panel defined by the layout function).
plot(NULL, xlim=c(start,end), ylim = c(0,200), axes=FALSE, ann=FALSE)
# Tell R the next plot call will be to the one that is already there and not initiate a new one.
par(new = TRUE)
# Plot the ATAC-seq track.
plot(0.5*(start(melp_5th_brain) + end(melp_5th_brain))+plotDiff, melp_5th_brain$score, type='l', xlim = c(start,end), ylim = c(0,200), ylab = "", yaxt = "n", lwd = 1, xlab = "", xaxt = "n", main = "", bty='none', col = "black")

# Add some text specifying the stage and tissue of the sample.
mtext('ATAC-seq 5th instar brain', side = 1, cex=0.8, padj = 0, las = 1, adj=1)
# Add the y-axis labels
axis(2, at = seq(0,200, by=50), line = 1)
mtext('Score', side = 2, cex=0.8, line = 3)

# Sample 2

plot(NULL, xlim=c(start,end), ylim = c(0,200), axes=FALSE, ann=FALSE)
par(new = TRUE)
plot(0.5*(start(melp_5th_FW) + end(melp_5th_FW))+plotDiff, melp_5th_FW$score, type='l', xlim = c(start,end), ylim = c(0,200), ylab = "", yaxt = "n", lwd = 1, xlab = "", xaxt = "n", main = "", bty='none', col = "black")

mtext('ATAC-seq 5th instar FW', side = 1, cex=0.8, padj = 0, las = 1, adj=1)
axis(2, at = seq(0,200, by=50), line = 1)
mtext('Score', side = 2, cex=0.8, line = 3)

# Sample 3

plot(NULL, xlim=c(start,end), ylim = c(0,200), axes=FALSE, ann=FALSE)
par(new = TRUE)
plot(0.5*(start(erato_5th_brain) + end(erato_5th_brain)), erato_5th_brain$score, type='l', xlim = c(start,end), ylim = c(0,200), ylab = "", yaxt = "n", lwd = 1, xlab = "", xaxt = "n", main = "", bty='none', col = "black")

mtext('ATAC-seq 5th instar brain', side = 1, cex=0.8, padj = 0, las = 1, adj=1)
axis(2, at = seq(0,200, by=50), line = 1)
mtext('Score', side = 2, cex=0.8, line = 3)

# Sample 4

plot(NULL, xlim=c(start,end), ylim = c(0,200), axes=FALSE, ann=FALSE)
par(new = TRUE)
plot(0.5*(start(erato_5th_FW) + end(erato_5th_FW)), erato_5th_FW$score, type='l', xlim = c(start,end), ylim = c(0,200), ylab = "", yaxt = "n", lwd = 1, xlab = "", xaxt = "n", main = "", bty='none', col = "black")

mtext('ATAC-seq 5th instar FW', side = 1, cex=0.8, padj = 0, las = 1, adj=1)
axis(2, at = seq(0,200, by=50), line = 1)
mtext('Score', side = 2, cex=0.8, line = 3)