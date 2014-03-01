#!/usr/bin/env Rscript

plot.rel.dist = function(path) {
  x = read.table(path, header=T, sep="\t")
  y = split(x, x$file)
  for(f in names(y)) {
    z = y[[f]]
    counts = table(z$centerdist)
    cnt = as.vector(counts)
    d = as.numeric(names(counts))
    print(d)
    print(cnt)
    plot(d, cnt, ylim=c(0, max(cnt)), xlab="Distance from center [nt]", ylab="Count", type='h', main=f)
  }
}

args = commandArgs(T)
print(args)
for(arg in args)
  plot.rel.dist(arg)
