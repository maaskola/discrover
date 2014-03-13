
ARGV.each{|path|
  out = path.gsub(/viterbi.*/, "posterior.txt")
  if out == path
    puts "Warning: file #{path} does not seem to be a Viterbi output file. Skipping."
    next
  end
  if File.exists?(out)
    puts "Warning: output file #{out} exists."
    next
  end
  puts "Writing to #{out}"
  o = File.open(out, "w")
  cmd = "echo #{path}"
  cmd = "zcat #{path}" if path =~ /.gz$/
  f = ""
  first = true
  `#{cmd}`.each_line{|line|
    if line =~ /^# (.*) details following/
      f = $1
                  # V-sites = 0/0 E-sites = 7.03064e-14/5.84726e-10 P(#sites>=1) = 7.81597e-14/5.8472e-10 Viterbi log-p = -60.7664
    elsif line =~ /^V-sites = (.*) E-sites = (.*) P\(#sites>=1\) = (.*) Viterbi log-p = (.*)/
      posterior = $3.split("/")
      if first
        o.print "path\t"
        posterior.each_with_index{|x,idx| o.print "\tV#{idx}" }
        o.puts
        first = false
      end
      o.puts f + "\t" + posterior.join("\t")
    end
  }
  `gzip -f #{out}`
}
