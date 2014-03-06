#!/usr/bin/env Rscript

plot.rel.dist = function(path) {
  x = read.table(path, header=T, sep="\t")
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

args = commandArgs(T)
print(args)
for(arg in args)
  plot.rel.dist(arg)
