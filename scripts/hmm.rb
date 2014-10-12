

def boostmatrix2matrix(x)
  matrix = []
  if x =~ /\[(\d+),(\d+)\](.*)/
    m = $1.to_i
    n = $2.to_i
    rest = $3
    rest.gsub(/^\(/,"").gsub(/\)$/,"").split("),(").map{|x|
      z = x.gsub(/^\(/,"").gsub(/\)$/,"").split(",").map{|y|
	y.to_f
      }
      matrix << z
    }
  end
  matrix
end

def boostify(x)
  n = x.length
  m = x[0].length
  y = "[#{n},#{m}]"
  y += "("
  y += x.map{|z|
    "(" + z.join(",") + ")"
  }.join(",")
  y += ")"
  y
end

$present2iupac = {
  "a" => "a",
  "c" => "c",
  "g" => "g",
  "t" => "t",
  "u" => "u",
  "ac" => "m",
  "gt" => "k",
  "ag" => "r",
  "ct" => "y",
  "at" => "w",
  "cg" => "s",
  "acg" => "v",
  "act" => "h",
  "agt" => "d",
  "cgt" => "b",
  "acgt" => "n"
}

$iupac2present = {
  "a" => "a",
  "c" => "c",
  "g" => "g",
  "t" => "t",
  "u" => "u",
  "m" => "ac",
  "k" => "gt",
  "r" => "ag",
  "y" => "ct",
  "w" => "at",
  "s" => "cg",
  "v" => "acg",
  "h" => "act",
  "d" => "agt",
  "b" => "cgt",
  "n" => "acgt"
}

def consensus_pos(x, threshold=0.1)
  present = []
  "acgt".split("").each_with_index{|y, idx| present << y if x[idx] > threshold}
  $present2iupac[present.join("")]
end

def consensus(x, threshold=0.1)
  x.map{|y| consensus_pos(y, threshold)}.join("")
end

module HMM
  class MotifStates
    attr_accessor :name, :states
    def initialize(name, states)
      @name = name
      @states = states
    end
  end
  class Parameters
    attr_accessor :transition, :emission, :first_state, :n_states, :n_emissions, :alpha, :beta, :version, :command, :time, :working_dir, :hmm_version, :plasma_version, :motifs, :classes, :order
    def initialize(f, verbose=false)
      @verbose = verbose
      line = f.gets.strip
      if line =~ /# HMM parameter format version (\d+)/
        @version = $1.to_i
        $stderr.puts "Version = #{version}" if verbose
        if @version == 1 or @version == 2 or @version == 3 or @version == 4 or @version == 5 or @version == 6
          line = f.gets.strip
          while line =~ /^#/
            if line =~ /# Command = (.*)/
              @command = $1
            elsif line =~ /# Run on (.*)/
              @time = $1
            elsif line =~ /# Run in (.*)/
              @working_dir = $1
            elsif line =~ /# Using plasma (.*)/
              @plasma_version = $1
            elsif line =~ /# hmm (.*)/
              @hmm_version = $1
            end
            line = f.gets.strip
          end
          # 4.times{|i| f.gets } if @version == 3
          case @version
          when 1 then
            @alpha = 1
            @beta = 0
          when 2..3 then
            if line =~ /Alpha (.*)/
              @alpha = $1.to_f
            else
              raise "Expecting 'Alpha ...'. Format conflict in line: #{line}"
            end
            line = f.gets.strip
            if line =~ /Beta (.*)/
              @beta = $1.to_f
            else
              raise "Format conflict in line: #{line}"
            end
            line = f.gets.strip
          end
          while line =~ /_prior/ or line =~ /^Dataset/
            line = f.gets.strip
          end
          if line =~ /(\d+) states/
            @n_states = $1.to_i
          else
            raise "Format conflict in line: #{line}"
          end
          $stderr.puts "n_states = #{@n_states}" if verbose

          line = f.gets.strip
          if line =~ /(\d+) emissions/
            @n_emissions = $1.to_i
          else
            raise "Format conflict in line: #{line}"
          end
          $stderr.puts "n_emissions = #{@n_emissions}" if verbose

          line = f.gets.strip

          case @version
          when 1..3 then
            if line =~ /Motif-start (\d+)/
              @first_state = $1.to_i
            else
              raise "Format conflict in line: #{line}"
            end
            $stderr.puts "first_state = #{@first_state}" if verbose

            if @version == 3
              line = f.gets.strip
              if line =~ /Motif-end (\d+)/
                @last_state = $1.to_i
              else
                raise "Format conflict in line: #{line}"
              end
              $stderr.puts "last_state = #{@last_state}" if verbose
            end

            @motifs = []
            @motifs << MotifStates.new("all", (first_state..last_state).to_a)

          when 4..6 then
            @motifs = []
            while line =~ /Motif "(.*)" (.*)/
              name = $1
              states = $2.split(" ").map{|i| i.to_i }
              @motifs << MotifStates.new(name, states)
              line = f.gets.strip
            end
            if line =~ /State class (.*)/
              @classes = $1.split(" ").map{|x| x.to_i}
            else
              raise "Format conflict in line: #{line}"
            end
            line = f.gets.strip
            if line =~ /State order (.*)/
              @order = $1.split(" ").map{|x| x.to_i}
            else
              raise "Format conflict in line: #{line}" unless line =~ /Transition matrix/
            end
          end

          line = f.gets.strip unless line =~ /Transition matrix/
          raise "Format conflict in line: #{line}" if line != "Transition matrix"
          @transition = []
          @n_states.times{|i|
            line = f.gets.strip
            j, rest = line.split(" ",2)
            j = j.to_i
            values = rest.split.map{|x| x.to_f}
            raise "Format conflict in line: #{line}" if i != j or values.length != @n_states
            @transition << values
          }

          line = f.gets.strip
          raise "Format conflict in line: #{line}" if line != "Emission matrix"

          @emission = []
          @n_states.times{|i|
            line = f.gets.strip
            j, rest = line.split(" ",2)
            j = j.to_i
            values = rest.split.map{|x| x.to_f}
            raise "Format conflict in line: #{line}" if i != j or values.length != @n_emissions
            @emission << values
          }

        else
          puts "Version \"#{@version}\" not known!"
        end
      else
        @first_state = first_line.to_i
        transition_s = f.gets
        emission_s = f.gets
        @transition = boostmatrix2matrix(transition_s)
        @emission = boostmatrix2matrix(emission_s)
        @n_states = @transition.size()
        @n_emissions = @emission[0].size()
      end
    end

    def to_s
      y = @first_state.to_s + "\n"
      y += boostify(@transition) + "\n"
      y += boostify(@emission) + "\n"
      y
    end

    def information_content(q, motif=nil)
      ic = 0
      case @version
      when 1..3 then
        (@first_state ... (@n_states)).to_a.each{|i|
          4.times{|j|
            p = @emission[i][j]
            ic += p * (Math::log(p) - Math::log(q[j])) if(p > 0)
          }
        }
      when 4 then
        idx = 0
        while motif.nil?
          motif = motifs.keys[idx] unless motifs.keys[idx] == "Constitutive"
          idx += 1
        end

        puts "Error: IC calculation is outdated."
        exit(-1)
        @motifs[motif].each{|i|
          4.times{|j|
            p = @emission[i][j]
            ic += p * (Math::log(p) - Math::log(q[j])) if(p > 0)
          }
        }
      else
        puts "Error: IC calculation not implemented for this HMM parameter file format."
        exit(-1)
      end
      ic = ic / Math::log(2.0)
      ic.round(3)
    end

    # Find those states reachable from the set of states given as argument.
    # NOTE: There's a slight hack in that transitions from the set of final
    # states to the set of initial states are not considered.
    def reachable_states(states, initial, final)
      # puts "calling reachable_states for #{states.join(",")}"
      r = {}
      states.each{|i|
        r[i] = []
        states.each{|j|
          r[i] << j if @transition[i][j] > 0 and not(final.include?(i) and initial.include?(j))# and j != states[0]
        }
      }
      # puts "in reachable_states -> #{r.map{|x,y| "#{x}:[#{y.join(",")}]"}.join(",")}"
      states.length.times{|n|
        r.each{|i,x|
          y = []
          x.each{|z|
            y << z
            r[z].each{|a|
              y << a
            }
            y = y.sort.uniq
          }
          r[i] = y
        }
      }
      # puts "called reachable_states -> #{r.map{|x,y| "#{x}:[#{y.join(",")}]"}.join(",")}"
      r
    end

    # Determine a topological order of the states of a motif.
    # If no motif is specified, use the first one.
    # Ignore transitions from set of final states to the set of initial ones.
    def topological_order(motif)
      if motif.nil?
        puts "Error: topological order needs a motif."
        exit(-1)
      end
      # motif = @motifs.keys[0] if motif.nil?
      states = motif.states
      initial = initial_states(motif)
      final = final_states(motif)
      reach = reachable_states(states, initial, final)
      states.sort{|x,y|
        if x == y
          0
        elsif reach[x].include?(y)
          -1
        else
          1
        end
      }
      # $stderr.puts "topo order = #{states}"
      states
    end

    def reachable_from_state(state)
      r = []
      @n_states.times{|i|
        r << i if @transition[state][i] > 0
      }
      r
    end

    def reachable_from_states(states)
      r = []
      states.each{|state|
        r += reachable_from_state(state)
      }
      r.sort.uniq
    end

    def initial_states(motif)
      bg = @motifs.find{|x| x.name == "Background" }
      reachable_from_bg = reachable_from_states(bg.states)
      initial = []
      motif.states.each{|state|
        initial << state if reachable_from_bg.include?(state)
      }
      # $stderr.puts "initial states of #{motif}: #{motif.name} = #{initial.to_s}"
      initial
    end

    def final_states(motif)
      final = []
      bg = @motifs.find{|x| x.name == "Background" }
      motif.states.each{|state|
        included = false
        bg.states.each{|bg_state|
          if @transition[state][bg_state] > 0
            included = true
            break
          end
        }
        final << state if included
      }
      # puts "final states of #{motif} = #{final.join(",")}"
      final
    end
  end

  class Motif
    attr_accessor :name, :iupac, :ic
    def initialize(name, iupac, ic)
      @name = name
      @iupac = iupac
      @ic = ic
    end
    def to_s
      "name = #{@name}, iupac = #{@iupac}, ic = #{@ic}"
    end
  end

  class Summary
    attr_accessor :path, :motifs, :scores
    def initialize(path)
      @path = path
      @motifs = []
      @scores = {}
      File.open(path){|f|
        f.gets # Motif.summary
        f.gets # Motif name  Consensus     IC [bit]
        while (l = f.gets.strip) != ""
          name, iupac, ic = l.split()
          ic = ic.to_f
          @motifs << Motif.new(name, iupac, ic)
        end
        f.each{|line|
          line.strip!
          processed = false
          @motifs.each{|motif|
            if line =~ /(.*) decoded (.*) counts - #{motif.name}/
              seqs, scores = parse_scores(motif.name, f)
              posterior_or_viterbi = $1
              motif_or_site = $2
              kind = (posterior_or_viterbi + " " + motif_or_site).to_sym
              @scores[motif.name] = {} unless @scores.has_key?(motif.name)
              @scores[motif.name][{:contrast => seqs, :kind => kind}] = scores
              processed = true
              break
            end
          }
          # $stderr.puts "Ignored: " + line unless processed
        }
      }
      @scores.each{|motif,scores|
        scores.each{|score|
        $stderr.puts [motif, score].join(" ")
        }
      }
    end
    def parse_scores(motif, f)
      scores = {}
      seqs = {}
      line = f.gets
      line = f.gets.strip
      until line =~ /.* - #{motif} /
        seq, present, absent, percent = line.split()
        seqs[seq] = {:present => present.to_f, :absent => absent.to_f, :percent => percent.to_f}
        line = f.gets
      end
      final = false
      while line =~ /.* - #{motif} /
        case(line)
        when /Discriminatory mutual information = (.*)/
          scores[:mi] = $1.to_f
        when /Expected discriminatory mutual information = (.*)/
          scores[:exp_mi] = $1.to_f
        when /Variance of discriminatory mutual information = (.*)/
          scores[:var_mi] = $1.to_f
        when /Std\. dev\. of discriminatory mutual information = (.*)/
          scores[:sd_mi] = $1.to_f
        when /Bonferroni corrected log-P\(Chi-Square\(G-Test\)\) = (.*)/
          scores[:bonferroni] = $1.to_f
          final = true
        when /Log-P\(Chi-Square\(G-Test\)\) = (.*)/
          scores[:gtest_pvalue] = $1.to_f
        when /P\(Chi-Square\(G-Test\)\) = (.*)/
          scores[:gtest_pvalue] = $1.to_f
        when /G-test = (.*)/
          scores[:gtest_statistic] = $1.to_f
        when /Matthews correlation coefficient = (.*)/
          scores[:mcc] = $1.to_f
        when /DIPS t-score = (.*)/
          scores[:dips_tscore] = $1.to_f
        when /DIPS site-score = (.*)/
          scores[:dips_sitescore] = $1.to_f
        when /log P correct classification = (.*)/
          scores[:correctclass_logp] = $1.to_f
        when /og P likelihood difference = (.*)/
          scores[:dlogl] = $1.to_f
        end
        break if final
        line = f.gets
        break if line.nil?
        line.strip!
      end
      [seqs, scores]
    end
  end
end

