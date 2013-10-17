#!/usr/bin/env ruby
# A program to extract n-mer statistics and generate according random sequences
require 'optparse'
require 'ostruct'
require 'bio'

class MyOptParser
  def self.parse(args)
    options = OpenStruct.new
    options.verbose = false
    options.order = 1
    options.seed = nil
    options.samples = 1
    opts = OptionParser.new {|opts|
      opts.banner = "Usage: markov1shuffle.rb [options] path ...\nThis program will estimate a first-order Markov chain for each sequence in the FASTA files, and will generate sequences of identical length to the input sequences from these Markov chains."
      opts.separator ""
      opts.on("-v","--[no-]verbose", "Run verbosely") {|x|
        options.verbose = x
      }
      opts.on("-n","--samples N",Integer,"Number of samples to create for each entry sequence. Default = 1") {|n|
        options.samples = n
      }
      opts.on("-s","--seed N",Integer,"Seed to initialize the pseudo random number generator.") {|k|
        options.seed = k
      }
    }
    opts.parse!(args)

    if options.order <= 0
      puts "Error: Order of the Markov chain must be non-negative."
      exit(-1)
    end

    if options.seed.nil?
      options.seed = Random.new_seed
    end

    $stderr.puts "Seed for pseudo random number generation = #{options.seed}"
    srand(options.seed)

    options
  end
end



def normalize(x)
  z = x.values.reduce(:+)
  return Hash[x.map {|key,val| [key,val/z] } ]
end

def random(x)
  z = rand()
  c = 0.0
  last_k = nil
  x.each {|key,val|
    c += val
    if c >= z
      return key
    end
    last_k = key
  }
  return last_k
end

class MarkovChain
  def initialize(seq, k, alphabet="acgt".split(""), pseudo_count = 0.001)
    @pseudo_count = pseudo_count
    @k = k
    @alphabet = alphabet
    @mono_nucl = Hash.new(0.0)
    @alphabet.each {|a| @mono_nucl[a] = @pseudo_count}
    default = Hash.new(0.0)
    @alphabet.each {|a| default[a] = 1.0}
    default = normalize(default)
    @p = Hash.new(default)
    (seq.length - @k).times {|i|
      s = seq[i...(i+@k)]
      t = seq[i+@k]
      @mono_nucl[t] += 1.0
      unless @p.has_key?(s)
        @p[s] = Hash.new(0)
        @p[s] = Hash[@alphabet.map {|a| [a,@pseudo_count]}]
      end
      # puts [i, t].join("\t")
      @p[s][t] += 1.0
    }
    @p.each {|from,d|
      @p[from] = normalize(d)
    }
    @mono_nucl = normalize(@mono_nucl)
  end
  def sample(n)
    # puts @p
    x = ''
    (@k).times {|i|
      x += random(@mono_nucl)
    }
    z = x
    (n-@k).times {|i|
      # puts x
      y = random(@p[x])
      z += y
      x = x[1...@k] + y
    }
    z
  end
end

options = MyOptParser.parse(ARGV)

ARGV.each {|path|
  Bio::FastaFormat.open(path).each {|entry|
    # puts entry
    mchain = MarkovChain.new(entry.seq.downcase, options.order, entry.seq.downcase.split("").sort.uniq)
    options.samples.times {|i|
      head = ">Shuffle"
      head += "_#{i}" if options.samples > 1
      puts "#{head}_#{entry.definition}"
      puts mchain.sample(entry.seq.length)
    }
  }
}

