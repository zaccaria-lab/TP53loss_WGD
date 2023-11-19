library(ggplot2)
library(mclust)
library(biomaRt)
library(reshape2)
library(pheatmap)
library(gridExtra)
library(grid)
library(DNAcopy)
library(xlsx)
library(RColorBrewer)
library(beeswarm)
library(ggbeeswarm)

outputdir <- "/camp/lab/swantonc/working/albakim/PC9.CellLine/CNVSimone/"
cellline.root <- list.files(outputdir, pattern = "PC9", full.names = F)
cellline.bin.dirs <- list.files(paste0(outputdir, cellline.root), full.names = T, recursive = F)
BAF <- 0.5
pval <- 0.05
min.total.depth <- 40
min.var.counts <- 10
resultsdir <- "/camp/lab/swantonc/working/albakim/PC9.CellLine/CNVSimone/results"
set.seed(7)
pdf.path <- "x"

pdf(pdf.path, width = 21, height = 7, onefile = TRUE)
summary <- NULL
for(i in 1:length(cellline.bin.dirs)){
        print(i)
        print(cellline.bin.dirs[i])
        derived.cl <- unlist(lapply(strsplit(cellline.bin.dirs[i], split = "/"), '[', 10))
        derived.cl <- gsub("\\.1Mb.*", "", derived.cl)
        progenitor.cl <- gsub("SU_.*", "SU_O0", derived.cl)
        
        ##read in the combined BIN/LogR df; colnames need correcting
        bin1M <- read.delim(cellline.bin.dirs[i], sep="\t", header = F, stringsAsFactors = F, skip=1)
        cols <- readLines(cellline.bin.dirs[i], n=1)
        cols <- gsub("pos\t", "pos\tpos.1\t", cols)
        cols <- unlist(strsplit(cols, split = "\t"))
        colnames(bin1M) <- cols
        ##now to make colnames not require the actual cell line names
        colnames(bin1M) <- gsub(derived.cl, "derived", colnames(bin1M))
        colnames(bin1M) <- gsub(progenitor.cl, "progenitor", colnames(bin1M))
        
        ##filter the file to remove SNPs with zero coverage
        bin1M <- bin1M[!bin1M$derived_tot%in%0,]
        bin1M <- bin1M[!bin1M$progenitor_tot%in%0,]
        
        ##filter the file to ensure minimum tumour SNP coverage as per defined thresholds
        bin1M <- bin1M[bin1M$derived_var>=min.var.counts & bin1M$derived_tot>=min.total.depth, ]
        
        ## remove the sex chromosomes
        bin1M <- bin1M[!bin1M$CHR%in%c("chrX", "chrY"),]
        
        ## now we need to define whether a SNP is significnatly different from BAF 0.5
        test <- function(x, n, p){binom.test(x, n, p, alternative="two.sided")}
        results <- mapply(test, bin1M$derived_var, bin1M$derived_tot, BAF)
        results <- t(results)
        results <- data.frame(results, stringsAsFactors = F)
        bin1M$pval <- unlist(results$p.value)
        bin1M$binom.pass.fail <- ifelse(bin1M$pval<pval, "FAIL", "PASS") ## PASS = at 0.5, FAIL != 0.5
        
        ## now we would like to get the depth ratio & corrected logR
        bin1M <- bin1M[!bin1M$REF_COUNTS%in%0 & !bin1M$TUMOUR_COUNTS%in%0, ]
        bin1M$read.depth.ratio <- bin1M$TUMOUR_COUNTS/bin1M$REF_COUNTS
        total.no.reads.tumour <- sum(bin1M$TUMOUR_COUNTS)
        total.no.reads.normal <- sum(bin1M$REF_COUNTS)
        correction.ratio <- total.no.reads.normal /total.no.reads.tumour 
        bin1M$read.depth.ratio.corr <- bin1M$TUMOUR_COUNTS/bin1M$REF_COUNTS *correction.ratio
        bin1M$LogR.corr <- log2(bin1M$read.depth.ratio.corr)
        
        #a <- ggplot() + geom_point(data = bin1M, aes(x=as.numeric(row.names(bin1M)), y= read.depth.ratio.corr, col=binom.pass.fail), alpha=0.2) +
        #        theme_minimal() + ggtitle(derived.cl)
        #b <- ggplot() + geom_point(data = bin1M, aes(x=as.numeric(row.names(bin1M)), y= derived_BAF, col=binom.pass.fail), alpha=0.2) +
        #        theme_minimal()
        #a.bar <- ggplotGrob(a + theme(plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm")))
        #b.bar <- ggplotGrob(b + theme(plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm")))
        #a.bar$widths <- unit.pmax(a.bar$widths, b.bar$widths)
        #b.bar$widths <- unit.pmax(a.bar$widths, b.bar$widths)
        
        #combined.list        <- list(a.bar, b.bar)
        #combined.plot.layout <- matrix(c(c(rep(1, 1)), c(rep(2, 1))), ncol = 1, byrow = T)
        #combined.plot        <- arrangeGrob(grobs = combined.list, nrow = nrow(combined.plot.layout), ncol = 1, layout_matrix = combined.plot.layout)
        
        #grid.draw(combined.plot)
        
        
        df <- bin1M[!duplicated(bin1M[c("chrom", "START")]),]
        cna <- CNA(genomdat=df$read.depth.ratio.corr, chrom=df$chrom, maploc=df$START, data.type='logratio')
        segs <- segment(cna, verbose=0)
        segs <- segs$output
        
        CNA.object <- CNA( genomdat = df$read.depth.ratio.corr, chrom = df$chrom, maploc = df$START, data.type = 'logratio')
        CNA.object.smoothed     <- smooth.CNA(CNA.object)
        CNA.object.smoothed.seg <- segment(CNA.object.smoothed, verbose=0, min.width=2)
        seg.pvalue <- segments.p(CNA.object.smoothed.seg, ngrid=100, tol=1e-6, alpha=0.05, search.range=100, nperm=1000)
        
        
        if(nrow(segs)<1){
                print("ERROR!!! not enough rows from DNAcopy")
        }
        
        segmeanbin1M <- rep(NA, nrow(bin1M))
        for(j in 1:nrow(segs)){
                rows.to.replace <- which(bin1M$chrom%in%segs$chrom[j] & bin1M$START>= segs$loc.start[j] & bin1M$START<= segs$loc.end[j])
                segmeanbin1M[rows.to.replace] <- segs$seg.mean[j]
        }
        bin1M$segmentedLogR <- segmeanbin1M
        
        #a <- ggplot() + geom_point(data = bin1M, aes(x=as.numeric(row.names(bin1M)), y= segmentedLogR, col=binom.pass.fail), alpha=0.2) +
        #        theme_minimal() + ggtitle(derived.cl)
        #b <- ggplot() + geom_point(data = bin1M, aes(x=as.numeric(row.names(bin1M)), y= derived_BAF, col=binom.pass.fail), alpha=0.2) +
        #        theme_minimal()
        #a.bar <- ggplotGrob(a + theme(plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm")))
        #b.bar <- ggplotGrob(b + theme(plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm")))
        #a.bar$widths <- unit.pmax(a.bar$widths, b.bar$widths)
        #b.bar$widths <- unit.pmax(a.bar$widths, b.bar$widths)
        
        #combined.list        <- list(a.bar, b.bar)
        #combined.plot.layout <- matrix(c(c(rep(1, 1)), c(rep(2, 1))), ncol = 1, byrow = T)
        #combined.plot        <- arrangeGrob(grobs = combined.list, nrow = nrow(combined.plot.layout), ncol = 1, layout_matrix = combined.plot.layout)
        
        #grid.draw(combined.plot)
        
        
        ## some bad bins may be noisy - therefore, ? only keep bins where >33% PASS
        bin1M$binID <- paste(bin1M$CHR, bin1M$START, bin1M$STOP, sep=":")
        uniquebins <- unique(bin1M$binID)
        usebin <- NULL
        meanBAFforbin <- NULL
        for(j in 1:length(uniquebins)){
                tmp <- bin1M[bin1M$binID%in%uniquebins[j],]
                perc.SNPs.pass <- length(which(tmp$binom.pass.fail%in%'PASS'))/nrow(tmp)*100
                t <- ifelse(perc.SNPs.pass>20, "use", "dontuse")
                pairs <- mapply(c, tmp$derived_var, tmp$derived_tot-tmp$derived_var, SIMPLIFY = FALSE)
                allBAFs <- unlist(lapply(pairs, min))/tmp$derived_tot
                to.keep <- which(allBAFs>0.2)
                if(length(to.keep)>=2){
                        meanBAFforbin <- c(meanBAFforbin, sum(unlist(lapply(pairs, min))[to.keep])/sum(tmp$derived_tot[to.keep]))
                }else{
                        meanBAFforbin <- c(meanBAFforbin,0)
                }
                usebin <- c(usebin, t)
        }
        bin1M$usebin <- ifelse(bin1M$binID%in%uniquebins[usebin%in%"use"], "use", "dontuse")
        table(bin1M$usebin)
        table(usebin)
        binBAF <- cbind(uniquebins, meanBAFforbin)
        binBAF <- data.frame(binBAF, stringsAsFactors = F)
        colnames(binBAF) <- c("binID", "meanBAFforbin")
        
        
        ## cluster the logR bins with SNPs at BAF 0.5 - the main/largest peak likely represents the logR0 andSNP0.5
        ## this means that this is the "balanced" cluster; i.e. no change from prog to derived
        bin1M.tmp <- bin1M[bin1M$usebin%in%"use",]
        #bin1M.tmp <- unique(bin1M.tmp[, c("read.depth.ratio.corr", "binID")])
        
        m <- Mclust(bin1M.tmp$segmentedLogR, modelNames = "E", G = 1:4)
        print(summary(m))
        #plot(m, what = "density")
        
        #pass.read.depth <-  cbind(bin1M.tmp$read.depth.ratio.corr, m$classification)
        pass.read.depth <-  cbind(bin1M.tmp, m$classification)
        pass.read.depth <- data.frame(pass.read.depth, stringsAsFactors = F)
        colnames(pass.read.depth)[ncol(pass.read.depth)]<- "cluster"
        
        #a <- ggplot(data = pass.read.depth, aes(x=as.numeric(row.names(pass.read.depth)), y= segmentedLogR, col=as.factor(cluster))) + geom_point( alpha=0.2) +
        #        theme_minimal() + ggtitle(derived.cl)
        #b <- ggplot(data = pass.read.depth, aes(x=as.numeric(row.names(pass.read.depth)), y= derived_BAF, col=as.factor(cluster))) + geom_point( alpha=0.2) +
        #        theme_minimal() + ggtitle(derived.cl)
        #a.bar <- ggplotGrob(a + theme(plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm")))
        #b.bar <- ggplotGrob(b + theme(plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm")))
        #a.bar$widths <- unit.pmax(a.bar$widths, b.bar$widths)
        #b.bar$widths <- unit.pmax(a.bar$widths, b.bar$widths)
        
        #combined.list        <- list(a.bar, b.bar)
        #combined.plot.layout <- matrix(c(c(rep(1, 1)), c(rep(2, 1))), ncol = 1, byrow = T)
        #combined.plot        <- arrangeGrob(grobs = combined.list, nrow = nrow(combined.plot.layout), ncol = 1, layout_matrix = combined.plot.layout)
        
        #grid.draw(combined.plot)
        
        
        summary.results <- summary(pass.read.depth$segmentedLogR)
        table(pass.read.depth$cluster)
        cluster.with.max.no.points <- as.numeric(names(table(pass.read.depth$cluster))[table(pass.read.depth$cluster)%in%max(table(pass.read.depth$cluster))])
        if(length(m$parameters$mean)==1){
                biggest.cluster.mean <- as.numeric(m$parameters$mean)
        }else{
                biggest.cluster.mean <- as.numeric(m$parameters$mean[names(m$parameters$mean)%in%cluster.with.max.no.points])
        }
        
        variance <- m$parameters$variance$sigmasq
        standard.dev <- sqrt(variance)
        
        
        #all.big.cluster.points.readdepth <- pass.read.depth$read.depth.ratio.corr[pass.read.depth$cluster%in%cluster.with.max.no.points]
        all.big.cluster.points.segmentedLogR <- pass.read.depth$segmentedLogR[pass.read.depth$cluster%in%cluster.with.max.no.points]
        left.ci <- quantile(all.big.cluster.points.segmentedLogR, probs = c(0.01, 0.99))[1]
        right.ci <- quantile(all.big.cluster.points.segmentedLogR, probs = c(0.01, 0.99))[2]
        
        #biggest.cluster.mean <- mean(unique(all.big.cluster.points.segmentedLogR))
        #standard.dev <- sd(unique(all.big.cluster.points.segmentedLogR))
        #n <- length(unique(all.big.cluster.points.segmentedLogR))
        #margin.error <- qt(0.95, n-1)*(standard.dev/sqrt(n))
        #left.ci <- biggest.cluster.mean - margin.error
        #right.ci <- biggest.cluster.mean + margin.error
        #hist(all.big.cluster.points.segmentedLogR)
        #abline(v = c(left.ci, right.ci), col="blue")
        
        
        z.test = function(a, mu, var){
                zeta = (mean(a) - mu) / (sqrt(var / length(a)))
                return(zeta)}
        
        #no.points <-  max(table(pass.read.depth$cluster))
        #z.test = function(a, mu, var){
        #        zeta = (mean(a) - mu) / (sqrt(var / no.points))
        #        return(zeta)}
        
        ReadD.zscore<-  lapply(
                bin1M$segmentedLogR ,
                FUN = function(x) {
                        z.test(x, biggest.cluster.mean, variance)
                })
        ReadD.zscore <- unlist(ReadD.zscore)
        ReadD.pvalue <- 2*pnorm(q = abs(ReadD.zscore),lower.tail = FALSE) ## two sided p value test
        bin1M$read.zscore <- ReadD.zscore
        bin1M$read.pvalue <- ReadD.pvalue
        bin1M$read.pvalue.right <- pnorm(q = ReadD.zscore,lower.tail = FALSE)
        bin1M$read.pvalue.left <- pnorm(q = ReadD.zscore,lower.tail = TRUE)
        
        
        ## define the sig gained and lost 
        bin1M$readD.sig <- ifelse(bin1M$read.pvalue<0.05, "sig", "not.sig")
        bin1M$gain.vs.loss <- ifelse(bin1M$segmentedLogR<biggest.cluster.mean, "loss", "gain")
        bin1M$final.status <- paste(bin1M$readD.sig, bin1M$gain.vs.loss, sep=":")
        bin1M$final.status <- ifelse(grepl("not", bin1M$final.status), "not.sig", bin1M$final.status)
        bin1M <- merge(bin1M, binBAF, by="binID", all.x=T, all.y=FALSE)
        
        bin1M$read.pvalue.group <- ifelse(bin1M$read.pvalue>0.05, "not.sig", NA)
        bin1M$read.pvalue.group <- ifelse(bin1M$read.pvalue>=0.001 & is.na(bin1M$read.pvalue.group), "slight.sig", bin1M$read.pvalue.group)
        bin1M$read.pvalue.group <- ifelse(bin1M$read.pvalue<0.001, "sig", bin1M$read.pvalue.group)
        
        #bin1M$read.pvalue.group <- ifelse(bin1M$read.depth.ratio.corr>right.ci | bin1M$read.depth.ratio.corr<left.ci, "sig", "not.sig")
        
        
        bin1M$meanBAFforbin <- as.numeric(bin1M$meanBAFforbin)
        
        write.table(bin1M, file = paste0(resultsdir, "/", derived.cl, "readdepth.txt"), sep="\t", col.names = T, row.names = F, quote = F)
        
        bin1M$CHR <- gsub("chr", "", bin1M$CHR)
        bin1M$CHR <- as.numeric(bin1M$CHR)
        bin1M$START <- as.numeric(bin1M$START)
        bin1M$STOP <- as.numeric(bin1M$STOP)
        bin1M <- bin1M[order(bin1M$CHR, bin1M$STOP),]
        rownames(bin1M) <- 1:nrow(bin1M)
        
        
        vlines <- NULL
        chr <- unique(bin1M$CHR)
        for(k in 1:length(chr)){
                vlines <- c(vlines, max(as.numeric(row.names(bin1M[bin1M$CHR%in%chr[k],]))))
        }
        
        labels <- vlines[1]/2
        for(k in 2: length(chr)){
                labels <- c(labels, ((vlines[k]-vlines[k-1])/2 + vlines[k-1]))
        }
        
        
        #ggplot(data = bin1M, aes(x=as.numeric(row.names(bin1M)), y= LogR.corr, col)) + geom_point(aes(col=cut(1-logR.pvalue, c(0, 0.95, 0.99, 1)))) +
        #theme_minimal() + scale_color_manual(values = c("light grey","yellow","red","black")) + ggtitle(derived.cl)
        a<- ggplot(data = bin1M, aes(x=as.numeric(row.names(bin1M)), y= segmentedLogR, col=read.pvalue.group)) + geom_point() +
                theme_classic() + scale_color_manual(values = c("light grey","red", "yellow")) + ggtitle(derived.cl) +
                theme() + xlab("position") + geom_vline(xintercept = vlines, col="light grey", lty=2) +
                scale_x_continuous(breaks = labels, labels = chr ) + ylim(c(0,3))
        #b <- ggplot(data = bin1M, aes(x=as.numeric(row.names(bin1M)), y= derived_BAF, col=read.pvalue.group)) + geom_point() +
        #        theme_classic() + scale_color_manual(values = c("light grey","red", "yellow")) + ggtitle(derived.cl) +
        #        theme() + xlab("position") + geom_vline(xintercept = vlines, col="light grey", lty=2) +
        #        scale_x_continuous(breaks = labels, labels = chr )
        #        
        #        
        #        #ggplot(data = bin1M, aes(x=as.numeric(row.names(bin1M)), y=meanBAFforbin, col=read.pvalue.group)) + geom_point() +
        #        #theme_minimal() + scale_color_manual(values = c("light grey","red","yellow")) + ggtitle(derived.cl)
        #a.bar <- ggplotGrob(a + theme(plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm")))
        #b.bar <- ggplotGrob(b + theme(plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm")))
        #a.bar$widths <- unit.pmax(a.bar$widths, b.bar$widths)
        #b.bar$widths <- unit.pmax(a.bar$widths, b.bar$widths)
        #
        #combined.list        <- list(a.bar, b.bar)
        #combined.plot.layout <- matrix(c(c(rep(1, 1)), c(rep(2, 1))), ncol = 1, byrow = T)
        #combined.plot        <- arrangeGrob(grobs = combined.list, nrow = nrow(combined.plot.layout), ncol = 1, layout_matrix = combined.plot.layout)
        
        print(a)
        
        
}
dev.off()
