
get_nucs <- function( sampd ){
  cat("__Processing chromosomes.\n")
  sampn <- unlist(lapply(strsplit(sampd, split = "//"), function(x){x[length(x)]}))
  #cmd <- paste("~/gits/fasta2nuccomp/fasta2nuccomp.py ", 
  #cmd <- paste("~/gits/nuccomp/nuccomp.py ", 
  cmd <- paste("~/gits/nuccomp/python/nuccomp.py ", 
               sampd,
               "/",
               sampn,
               #               "/",
               ".softmasked.fasta.gz",
               sep = "")
  
  if( !file.exists( paste(sampn, ".softmasked_nuccomp.csv", sep = "") ) ){
    cat( "Running nuccomp: ", sampn, ".\n" )
    system( cmd )
  }
  
  nucs <- read.csv( paste( sampn, ".softmasked_nuccomp.csv", sep = "") )
  nucs <- nucs[ grep( paste(sampn, ".chr", sep = ""), nucs$Id ), ]
  
  if( length( grep("chrY", nucs$Id) ) == 0 ){
    last_chrom <- nrow(nucs)
    nucs <- rbind(nucs, NA)
    nucs[ last_chrom + 1, -1] <- 0
    nucs[ last_chrom + 1, 1] <- "chrY"
    nucs[ last_chrom + 1, 2] <- 1
  }
  if( length( grep("chrX", nucs$Id)) == 0 ){
    last_chrom <- nrow(nucs)
    nucs <- rbind(nucs, nucs[last_chrom, ])
    nucs[ last_chrom + 0, -1] <- 0
    nucs[ last_chrom + 0, 1] <- "chrX"
    nucs[ last_chrom + 0, 2] <- 1
  }
  return( nucs )
}


##### ##### ##### ##### #####


get_wins <- function( sampd ){
  sampn <- unlist(lapply(strsplit(sampd, split = "//"), function(x){x[length(x)]}))
  
  cat("__Processing windows.\n")
  
  cmd <- paste("~/gits/hempy/bin/fast_win.py ", 
               sampd ,
               "/",
               sampn,
               #               "/",
               ".softmasked.fasta.gz",
               " --win_size 1000000",
               sep = "")
  
  if( !file.exists( paste(sampn, ".softmasked_wins.csv", sep = "") ) ){
    cat( "Running nuccomp: ", sampn, ".\n" )
    system( cmd )
  }
  
  wins <- read.csv( paste( sampn, ".softmasked_wins.csv", sep = "") )
  wins <- wins[ grep( paste(sampn, ".chr", sep = ""), wins$Id ), ]
  wins$chrn <- sub(".+chr", "", wins$Id)
  wins$chrn[ wins$chrn == "X" ] <- 10
  wins$chrn[ wins$chrn == "Y" ] <- 11
  wins$chrn <- as.numeric(wins$chrn)
  
  # Scale and center.
  # Scaline of CG windows (chromosome width) is on a per chromosome basis.
  #  wins$CGs <- 0
  wins$CGs <- wins$CG / wins$Win_length * 100
  #  wins$notCG <- (wins$Win_length / 1) - wins$CG
  #  wins$notCGs <- 0
  for( j in unique(wins$Id) ){
    wins$CGs[ wins$Id == j] <- wins$CGs[ wins$Id == j] - min(wins$CGs[ wins$Id == j], na.rm = TRUE)
    wins$CGs[ wins$Id == j] <- wins$CGs[ wins$Id == j] / max(wins$CGs[ wins$Id == j], na.rm = TRUE)
    #
    # wins$CGs[ wins$Id == j] <- wins$CG[ wins$Id == j] - min(wins$CG[ wins$Id == j], na.rm = TRUE)
    # wins$CGs[ wins$Id == j] <- wins$CGs[ wins$Id == j]/max(wins$CGs[ wins$Id == j], na.rm = TRUE)
    #
    # wins$notCGs[ wins$Id == j] <- wins$notCG[ wins$Id == j] - min(wins$notCG[ wins$Id == j], na.rm = TRUE)
    # wins$notCGs[ wins$Id == j] <- wins$notCGs[ wins$Id == j]/max(wins$notCGs[ wins$Id == j], na.rm = TRUE)
  }
  wins$iCGs <- 1 - wins$CGs
  
  #wins$ATs <- 1 - wins$CGs
  wins$ATs <- wins$Win_length/1e6 - wins$CGs
  
  #sampd <- "/media/knausb/E737-9B48/releases/scaffolded/AH3Ma/"
  #sampn <- "AH3Ma"
  #sampd <- "/media/knausb/E737-9B48/releases/scaffolded/GRMb/"
  #sampn <- "GRMb"
  #sampd <- "/media/knausb/E737-9B48/releases/scaffolded/SODLb/"
  #sampn <- "SODLb"
  genes <- read.table( 
    paste(sampd, "/", sampn, ".primary_high_confidence.gff3.gz", sep = "" ), 
    sep = "\t" )
  genes <- genes[genes[, 3] == "gene", ]
  # Select chromosomes, omit unplaced sequences.
  genes <- genes[grep("chr", genes[, 1]), ]
  #table(genes[, 1])
  
  # Windowize
  wins$gcnt <- 0
  for(i in 1:nrow(wins)){
    tmp <- genes[genes$V1 == wins$Id[i] & genes$V4 >= wins$Start[i] & genes$V5 < wins$End[i], ]
    wins$gcnt[i] <- nrow(tmp)
  }
  #hist(wins$gcnt)
  #
  wins[ wins$gcnt > 300, ]
  #hist(wins$gcntsc)
  
  # Scaling of gene count windows is on a per assembly basis.
  # wins$gcntsc <- wins$gcnt - min(wins$gcnt)
  # wins$gcntsc <- wins$gcnt - min(wins$gcnt)
  wins$gcntsc <- wins$gcnt/wins$Win_length
  wins$gcntsc <- wins$gcntsc - min(wins$gcntsc)
  wins$gcntsc <- wins$gcntsc / max(wins$gcntsc)
  my_index <- round( wins$gcntsc * 100 )
  my_index[ my_index <= 0] <- 1
  
  # Color ramp.
  #wins$gcol <- heat.colors(n=100)[ my_index ]
  #wins$gcol <- colorRampPalette(c("yellow", "orange", "red"))( 100 )[ my_index ]
  #wins$gcol <- colorRampPalette(c("red", "orange", "yellow"))( 100 )[ my_index ]
  #wins$gcol <- colorRampPalette(c("#0000FF", "#228B22", "#A0522D"))( 100 )[ my_index ]
  #wins$gcol <- colorRampPalette(c("#87CEEB", "#3CB371", "#228B22", "#A0522D"))( 100 )[ my_index ]
  
  #wins$gcol <- viridisLite::plasma(n = 100, alpha = 1, begin = 0.1, end = 1)[ my_index ]
  wins$gcol <- viridisLite::magma(n = 100, alpha = 1, begin = 0.2, end = 1.00)[ my_index ]
  
  return(wins)
}


##### ##### ##### ##### #####


plot_ideo <- function( sampd ) {
  suppressPackageStartupMessages(require(ggplot2))
  
  sampn <- unlist(lapply(strsplit(sampd, split = "//"), function(x){x[length(x)]}))
  cat( "Processing ", sampn, ".\n", sep = "" )
  nucs <- get_nucs( sampd )
  wins <- get_wins( sampd )
  
  blstn <- read.csv( 
    #    paste("/media/knausb/Vining_lab/knausb/blast_projects/todd_rrna/blastn/", sampn, "_blastn.csv", sep = ""),
    # paste("/media/knausb/Vining_lab/knausb/blast_projects/todd_rrna/blastn_uc/", sampn, "_blastn.csv", sep = ""),
    # paste("/media/knausb/E737-9B48/knausb/blast_projects/todd_rrna/blastn_uc/", sampn, "_blastn.csv", sep = ""),
    # header = FALSE)
    paste("../blastn_markers/", sampn, "_blastn.csv.gz", sep = ""),
    header = FALSE)
  colnames(blstn) <- c("qseqid","qlen","sseqid","slen","qstart","qend",
                       "sstart","send","evalue","bitscore","score","length",
                       "pident","nident","mismatch","positive","gapopen",
                       "gaps","ppos","sstrand") 
  #
  #length(grep("chr", blstn$sseqid, invert = TRUE))
  
  blstn <- blstn[grep("chr", blstn$sseqid), ]
  blstn$chrn <- sub(".*chr", "", blstn$sseqid)
  blstn$chrn[blstn$chrn == "X"] <- 10
  blstn$chrn[blstn$chrn == "Y"] <- 11
  blstn$chrn <- as.numeric(blstn$chrn)  
  
  # Orient chromosomes.
  # orient <- read.table("/media/knausb/E737-9B48/releases/scaffolded_csat_orientations.tsv", 
  #                      header = TRUE, sep = "\t")
  orient <- read.table("../ideograms_scaffolded/scaffolded_csat_orientations.tsv", 
             header = TRUE, sep = "\t")
  orient <- orient[orient$Sample == sampn, ]
  row.names(orient) <- paste(orient$Sample, orient$Chromosome, sep = ".")
  # orient[1:3, ]
  for(i in 1:nrow(orient)){
    if( orient$Flip[i] == "True" ){
      # nucs =  no action. Contains nucs$Length
      chrom_len <- nucs$Length[ nucs$Id == rownames(orient)[i] ]
      tmp <- wins[ wins$Id == rownames(orient)[i], ]
      #tmp <- wins[ wins$Id == rownames(orient)[i], , drop = FALSE]
      # tmp[1:3, 1:6]
      tmp$Start <- chrom_len - tmp$Start
      tmp$End <- chrom_len - tmp$End
      wins[ wins$Id == rownames(orient)[i], ] <- tmp
      #wins[ wins$Id == rownames(orient)[i], , drop = FALSE] <- tmp
      #blstn
      tmp <- blstn[ blstn$sseqid == rownames(orient)[i], , drop = FALSE]
      tmp[1:3, 1:8]
      tmp$sstart <- chrom_len - tmp$sstart
      tmp$send <- chrom_len - tmp$send
      blstn[ blstn$sseqid == rownames(orient)[i], ] <- tmp
    }
  }
  
  # blstn[1:3, ]
  # table(blstn$qseqid)
  
  chr_df <- data.frame(
    start = 1,
    end = nucs$Length,
    chr = nucs$Id
  )
  #  chr_df$chrn <- sub(".+chr", "", nucs$Id)
  chr_df$chrn <- sub(".*chr", "", nucs$Id)
  chr_df$chrn[chr_df$chrn == "X"] <- 10
  chr_df$chrn[chr_df$chrn == "Y"] <- 11
  chr_df$chrn <- as.numeric(chr_df$chrn)
  
  chrom_wid <- 0.02
  p <- ggplot2::ggplot()
  p <- p + ggplot2::geom_rect( 
    data = chr_df, 
    ggplot2::aes( xmin = chrn - chrom_wid,
                  xmax = chrn + chrom_wid,
                  ymin = start, ymax = end),
    fill = "#DCDCDC",
    color = "#000000"
  )
  
  thinw <- 0.28
  p <- p + ggplot2::geom_rect( 
    data = wins, 
    ggplot2::aes(
      #      xmin = chrn - ATs * thinw,
      #      xmax = chrn + ATs * thinw,
      xmin = chrn - iCGs * thinw,
      xmax = chrn + iCGs * thinw,
      # xmin = chrn - notCGs * thinw,
      # xmax = chrn + notCGs * thinw,
      ymin = Start, 
      ymax = End),
    fill = wins$gcol,
    color = NA
  )
  
  # Theme
  chr_names <- sub("^.+\\.", "", chr_df$chr)
  p <- p + ggplot2::scale_x_continuous( 
    breaks = chr_df$chrn,
    #labels = chr_df$chr
    labels = chr_names
  )
  
  p <- p + scale_y_continuous(
    breaks = seq( 0, 120e6, by = 10e6), 
    labels = seq( 0, 120, by = 10)
  )
  
  p <- p + ggplot2::theme_bw() + 
    ggplot2::theme( 
      #      panel.grid.minor.x = ggplot2::element_blank(), 
      panel.grid.minor.x = ggplot2::element_line( linewidth = 0.4, color = "#C0C0C0", linetype = 3 ),
      axis.text.x = element_text(angle = 60, hjust = 1),
      axis.title.x=element_blank(),
      panel.grid.major.y = ggplot2::element_line( linewidth = 0.4, color = "#C0C0C0", linetype = 1 ),
      panel.grid.minor.y = ggplot2::element_line( linewidth = 0.4, color = "#C0C0C0", linetype = 3 )
    )
  
  p <- p + ggplot2::ylab("Position (Mbp)")
  p <- p + ggtitle( sampn )
  
  my_locus <- ""
  table(blstn$qseqid)
  grep("5S_", blstn$qseqid)
  
  
  # Cannabinoids
  my_rows <- grep("AB2|LY|NC_044376.1:c427", blstn$qseqid)
  blstw <- 0.48
  p <- p + ggplot2::geom_rect( 
    data = blstn[my_rows, ], 
    ggplot2::aes( 
      xmin = chrn - blstw,
      xmax = chrn + blstw,
      ymin = sstart, 
      ymax = send),
    fill = "#228B22",
    color = '#228B22'
  )
  
  # 5S
  my_rows <- grep("5S_", blstn$qseqid)
  blstw <- 0.48
  p <- p + ggplot2::geom_rect( 
    data = blstn[my_rows, ], 
    ggplot2::aes( 
      xmin = chrn - blstw,
      xmax = chrn + blstw,
      ymin = sstart, 
      ymax = send),
    fill = "#000000",
    color = '#000000'
  )
  
  # 45S
  my_rows <- grep("45s_|26s_|5.8s_|18s_", blstn$qseqid)
  blstw <- 0.48
  p <- p + ggplot2::geom_rect( 
    data = blstn[my_rows, ], 
    ggplot2::aes( 
      xmin = chrn - blstw,
      xmax = chrn + blstw,
      ymin = sstart, 
      ymax = send),
    fill = "#B22222",
    linewidth = 2,
    #    linewidth = 4,
    color = '#B22222'
  )
  
  # # Cent, subteol
  # my_rows <- grep("AH3Ma_CS-237_satelite|CsatSD_centromere_237bp|CsatSD_centromere_370bp", blstn$qseqid)
  # #  blstw <- 0.48
  # p <- p + ggplot2::geom_rect( 
  #   data = blstn[my_rows, ], 
  #   ggplot2::aes( 
  #     xmin = chrn - 0.5,
  #     xmax = chrn - 0.2,
  #     ymin = sstart, 
  #     ymax = send),
  #   fill = "#1E90FF",
  #   color = '#1E90FF'
  # )
  
  # Cent, subteol
  my_rows <- grep("CsatSD_centromere_237bp", blstn$qseqid)
  #  blstw <- 0.48
  p <- p + ggplot2::geom_rect( 
    data = blstn[my_rows, ], 
    ggplot2::aes( 
      xmin = chrn - 0.5,
      xmax = chrn - 0.2,
      ymin = sstart, 
      ymax = send),
    fill = "#1E90FF",
    color = '#1E90FF'
  )
  
  # Cent, subteol
  my_rows <- grep("CsatSD_centromere_370bp", blstn$qseqid)
  #  blstw <- 0.48
  p <- p + ggplot2::geom_rect( 
    data = blstn[my_rows, ], 
    ggplot2::aes( 
#      xmin = chrn - 0.5,
#      xmax = chrn - 0.2,
      xmin = chrn + 0.2,
      xmax = chrn + 0.5,
      ymin = sstart, 
      ymax = send),
    fill = "#FF00FF",
    color = '#FF00FF'
  )
  
  return(p)
}

##### ##### ##### ##### #####





