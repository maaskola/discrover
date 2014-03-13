
init = function() {
  library(tikzDevice)
  library(caTools)
  library(png)
}

plot.ranks = function(x, k=500, width=3.25*1.0, height=3.5*1.0, path=NULL, dpi=1200) {
# plot.ranks = function(x, k=500, width=3.25, height=4, path=NULL, dpi=1200) {
  if(!is.null(path)) {
    tikz(path, width=width, height=height, standAlone=T, packages="\\usepackage{times,tikz}
\\usepackage[active,tightpage,psfixbb]{preview}
\\PreviewEnvironment{pgfpicture}
\\setlength\\PreviewBorder{0pt}
")
    # pdf(path, width=width, height=height)
  }
  y = split(x, x$path)
  n = length(y)
  lwd = 5
  cex = 0.8
  par(mfrow=c(n/2, 1))
  par(mar=c(3,4.5,1,0) + 0.1, cex=cex)
  for(i in 0:(n/2 - 1)) {
    if(dpi>0) {
      plot(runmean(y[[2*i+1]]$p_site,k), ylim=c(0,0.4), xlab="", ylab="", type='n')
      # mtext("bla", 2, 2.5)
      mtext("Motif presence", 2, 3.5, cex=cex)
      mtext(paste("averaged over",k,"sequences"), 2, 2.5, cex=cex)
      mtext(paste("Ranked sequences, RBM10 PAR-CLIP dataset", (i+1)), 1, 2, cex=cex)

      png.path = paste("sourceimg_", i+1, ".png", sep="")
      # tiff.path = paste("/tmp/blabla_", i, ".tiff", sep="")
      plt = par()$plt
      w = width * (plt[2] - plt[1]) # width of the plot region in inches
      h = height / (n/2) * (plt[4] - plt[3]) # height of the plot region in inches
      png(png.path, width=dpi*w, height=dpi*h)
      # tiff(tiff.path, width=dpi*w, height=dpi*h)
      par(mar=c(0,0,0,0), bty="n")
      plot(runmean(y[[2*i+1]]$p_site,k), type='l', ylim=c(0,0.4), ylab="", xlab="", lwd=lwd)
      lines(runmean(y[[2*i+2]]$p_site, k), col='red', lwd=lwd)
      dev.off()
      ima = readPNG(png.path)
      lim <- par()$usr
      rasterImage(ima, lim[1], lim[3], lim[2], lim[4])
      box()
    } else {
      plot(runmean(y[[2*i+1]]$p_site,k), type='l', ylim=c(0,0.4), ylab="",xlab="")
      lines(runmean(y[[2*i+2]]$p_site, k), col='red')
    }
  }
  if(!is.null(path))
    dev.off()
}


# function(path) {
#   #Replace the directory and file information with your info
#   ima <- readPNG(path)
# 
#   #Set up the plot area
#   plot(1:2, type='n', main="Plotting Over an Image", xlab="x", ylab="y")
# 
#   #Get the plot information so the image will fill the plot box, and draw it
#   lim <- par()
#   rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
# }
