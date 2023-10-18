##################################################################################
#### Written by T. Sauvage (2022). Ifremer, Centre Atlantique, Nantes, France ####
##################################################################################

trop <- function (datfile,seqlengths) {
# Import as list because of variable column number, and parse into usable data.table format
message("Importing and parsing trf's file")
dat <- scan(datfile, what = "", sep = "\n", quiet = TRUE) 
# PS: Switch the quiet argument above to false to see the number of lines read
pos <- regexpr("Sequence|^[[:digit:]]+\\s[[:digit:]]+",dat)
pos <- regmatches(dat, pos)
dna <- regexpr("Sequence:\\s.{36}|^[[:digit:]]+\\s[[:digit:]]+",dat)
dna <- regmatches(dat, dna)
rm(dat)
pos <- strsplit(pos, "Sequence")
pos <- lapply(pos, function(x){x[!x == ""]})
pos <- data.frame(do.call("rbind", strsplit(as.character(pos), " ")))
pos[pos == "character(0)"] <- "NA"
pos <- data.table(pos)
colnames(pos) <- c("startpos","endpos")
dna <- strsplit(dna, "Sequence: ")
dna <- lapply(dna, function(x){x[!x == ""]})
dna <- data.frame(do.call("rbind", strsplit(as.character(dna), " ")))
x <- which(grepl('^-?(0|[1-9][0-9]*)(\\.[0-9]+)?$',dna$X2))	# find numerical
dna$X2[x] <- NA	# replace by NAs
dna$X1 <- NULL
rm(x)
dna <- data.table(dna)
dna <- dna[, seqname := X2[1], .(cumsum(!is.na(X2)))]	# making new column seems important, fill in NA with sequence name
dna$X2 <- NULL
dna <- cbind(dna$seqname,pos)
rm(pos)
dna <- unique(dna)
dna$startpos <- gsub('NA', "0", dna$startpos)
dna$endpos <- gsub('NA', "0", dna$endpos)
dna <- dna[, startpos := as.numeric(startpos)]
dna <- dna[, endpos := as.numeric(endpos)]
colnames(dna) <- c("seqname","startpos","endpos")
# VERY necessary!!! order startpos increasing per sequence to make sure trgroup form properly below
dna <- dna[, setorder(.SD, startpos), by = seqname]
message("Found tandems from ", length(unique(dna$seqname))," sequences")
message("Merging overlapping tandems...")
# Create trgroup column according to interval overlap (i.e. different number show non-overlap) per read
dna[ , trgroup := cumsum(cummax(shift(endpos, fill = endpos[1])) < startpos) + 1, , by = seqname]
dna$trgroup <-(dna$trgroup)-1
dna[, minpos := min(startpos), by = list(seqname,trgroup)]
dna[, maxpos := max(endpos), by = list(seqname,trgroup)]
dna[,':='(startpos= NULL, endpos = NULL)]
dna <- unique(dna)
# Compute each repeat length
dna[, trlength := (maxpos-minpos)]
# Export
write.table(dna, "0_tandem_groups_per_reads.txt" , row.names = FALSE, quote = FALSE)
# Compute sum, max and merge
dna1 <- dna[, sum(trlength),by = seqname]
dna2 <- dna[, max(trgroup),by = seqname]
rm(dna)
dna3 <- merge(dna2, dna1, by.x = "seqname", by.y ="seqname")
rm(dna1)
rm(dna2)
setnames(dna3, c("seqname","trgroup","trlength"))
message("Importing sequence lengths")
seqlengths <- read.table(seqlengths, header = FALSE, sep="\t")
seqlengths <- as.data.table(seqlengths)
colnames(seqlengths) <- c("seqname","seqlength")
message("Found ", length(unique(seqlengths$seqname))," sequences")
seqlengths$seqlength <- as.numeric(seqlengths$seqlength)
message("Processing...")
dna3 <- merge(dna3,seqlengths, by.x = "seqname",by.y = "seqname", all.x = TRUE)
dna3$trprop <- round((dna3$trlength*100)/dna3$seqlength, digits = 2)
dna3 <- dna3[order(seqname),]
# Export tandem sum per read and keep number of fragments recorded per read 
# (i.e. = trgroup, which is non-overlapping tandem island/region per read, i.e or frequency)
write.table(dna3, "1_tandem_length_per_reads.txt" , row.names = FALSE, quote = FALSE)
message("File 1 exported")
# Compute summary stats
# 1 For number of sequences
total <- length(dna3$seqname)
no_tandem <- length(which(dna3$trlength == 0))
tandem <- total - no_tandem
mycounts <- rbind(no_tandem,tandem,total)
mycounts <- data.table(mycounts, keep.rownames = TRUE)
proptr <- round((tandem*100)/total, digits = 2)
propnotr <- round((no_tandem*100)/total, digits = 2)
propall <- round((total*100)/total, digits = 2)
myprops <- rbind(propnotr,proptr,propall)
myprops <- data.table(myprops)
mysummary <- cbind(mycounts,myprops$V1)
# 2 For cumulative base pair
total <- sum(dna3$seqlength)
tandem <- sum(dna3$trlength)
no_tandem <- total - tandem
mycum <- rbind(no_tandem,tandem,total)
mycum <- data.table(mycum, keep.rownames = TRUE)
proptr_cum <- round((tandem*100)/total, digits = 2)
propnotr_cum <- round((no_tandem*100)/total, digits = 2)
propall_cum <- round((total*100)/total, digits = 2)
myprops_cum <- rbind(propnotr_cum,proptr_cum,propall_cum)
myprops_cum <- data.table(myprops_cum)
mysummary <- cbind(mysummary,mycum$V1,myprops_cum$V1)
colnames(mysummary) <- c("class","nseq","nseq_perc","bp","bp_perc")
# Export summary
write.table(mysummary, "3_tandem_output_summary.txt" , row.names = FALSE, quote = FALSE)
# Get non 0 rows from tandem group per read, append seqlengths and compute pos
dna <- read.table("0_tandem_groups_per_reads.txt", header = TRUE)
dna <- data.table(dna)
dna <- dna[dna$trlength != 0, ]
dna <- merge(dna,seqlengths, by.x = "seqname",by.y = "seqname", all.x = TRUE)
rm(seqlengths)
dna$trprop <- round((dna$trlength*100)/dna$seqlength, digits = 2)
dna$start_3prime <- ((dna$seqlength-dna$maxpos)*(-1))-1
dna$stop_3prime <- ((dna$seqlength-dna$minpos)*(-1))-1
# Get 0 tandem rows only from summed file and add columns
dna3 <- dna3[dna3$trlength == 0, ]
dna3$minpos <- 0 
dna3$maxpos <- 0
dna3$start_3prime <- 0
dna3$stop_3prime <- 0
dna3 <- rbindlist(list(dna,dna3),use.names=TRUE)
rm(dna)
dna3 <- dna3[order(seqname),]
names(dna3)[names(dna3)=="minpos"] <- "start_5prime"
names(dna3)[names(dna3)=="maxpos"] <- "stop_5prime"
# Export and clean
write.table(dna3, "2_tandem_locations_per_read.txt" , row.names = FALSE, quote = FALSE)
message("File 2 exported")
message("File 3 exported")
rm(dna3)
unlink("0_tandem_groups_per_reads.txt")
}
