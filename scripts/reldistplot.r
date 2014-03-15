#!/usr/bin/env Rscript

plot.rel.dist = function(x) {
  y = split(x, x$file)
  for(f in names(y)) {
    z = y[[f]]
    m = split(z, z$motifidx)
    idx = 1
    for(g in names(m)) {
      n = m[[g]]
      counts = table(n$centerdist)
      cnt = as.vector(counts)
      d = as.numeric(names(counts))
      print(d)
      print(cnt)
      main = f
      if(length(m) > 1)
        main = paste(f, idx)
      idx = idx + 1
      plot(d, cnt, ylim=c(0, max(cnt)), xlab="Distance from center [nt]", ylab="Count", type='h', main=main)
    }
  }
}


plot.rel.dist.two = function(x, mdist=250, cex=2) {
  par(mar=c(4,3.5,0.5,0)+0.1, cex=cex)
  y = split(x, x$file)
  nfiles = length(y)
  if(length(y) != 2) {
    print("Error: this function assumes that three are just two files.")
    return()
  }
  signal = grep("signal", names(y))
  control = grep("signal", names(y), invert=T)
  names(y)[signal] = "signal"
  names(y)[control] = "control"

  signal = y[["signal"]]
  control = y[["control"]]

  m.signal = split(signal, signal$motifidx)
  m.control = split(control, control$motifidx)
  # idx = 1
  for(g in names(m.signal)) {
    n.signal = m.signal[[g]]
    n.control = m.control[[g]]

    counts.signal = table(n.signal$centerdist)
    counts.control = table(n.control$centerdist)

    cnt.signal = as.vector(counts.signal)
    cnt.control = as.vector(counts.control)

    mcnt = max(cnt.signal, cnt.control)

    d.signal = as.numeric(names(counts.signal))
    d.control = as.numeric(names(counts.control))

    md = max(mdist, abs(d.signal), abs(d.control))

    # print(d.signal)
    # print(cnt.signal)
    # print(d.control)
    # print(cnt.control)

    # main = f
    # if(length(m) > 1)
    #   main = paste(f, idx)
    main = ""
    # main = idx
    # idx = idx + 1
    d = c(d.signal, d.control)
    cnt = c(cnt.signal, cnt.control)
    col = c(rep("black", length(d.signal)), rep("red", length(d.control)))
    # plot(d, cnt, xlim=c(-md,md), ylim=c(0, mcnt), xlab="Distance from center [nt]", ylab="Count", type='h', main=main, col=col)
    plot(d, cnt, xlim=c(-md,md), ylim=c(0, mcnt), type='h', main=main, col=col, xlab="", ylab="")
    mtext("Distance from center [nt]", side=1, line=2.5, cex=cex)
    mtext("Count", side=2, line=2.5, cex=cex)
  }
}


rel.dist = function(path, out=NULL) {
  x = read.table(path, header=T, sep="\t")
  if(!is.null(out))
    pdf(out, width=10, height=5)
  plot.rel.dist.two(x, cex=2)
  if(!is.null(out))
    dev.off()
}

args = commandArgs(T)
# print(args)
for(arg in args) {
  print(arg)
  out = gsub(".table.*", ".pdf", arg)
  if(arg == out)
    print("Error: output file naming failed.")
  else
    rel.dist(arg, out=out)
}
