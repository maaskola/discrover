

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

module HMM
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

            @motifs = {}
            @motifs["all"] = (first_state..last_state).to_a

          when 4..6 then
            @motifs = {}
            while line =~ /Motif "(.*)" (.*)/
              name = $1
              states = $2.split(" ").map{|i| i.to_i }
              @motifs[name] = states
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
              raise "Format conflict in line: #{line}"
            end
          end

          line = f.gets.strip
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
        w = @n_states - @first_state
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

        @motifs[motif].each{|i|
          4.times{|j|
            p = @emission[i][j]
            ic += p * (Math::log(p) - Math::log(q[j])) if(p > 0)
          }
        }
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
    def topological_order(motif=nil)
      motif = @motifs.keys[0] if motif.nil?
      states = @motifs[motif]
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
      reachable_from_bg = reachable_from_states(@motifs["Background"])
      initial = []
      @motifs[motif].each{|state|
        initial << state if reachable_from_bg.include?(state)
      }
      # puts "initial states of #{motif} = #{initial.join(",")}"
      initial
    end

    def final_states(motif)
      final = []
      @motifs[motif].each{|state|
        included = false
        @motifs["Background"].each{|bg_state|
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

  class Summary
    attr_accessor :path, :mi, :gtest_statistic, :gtest_pvalue, :gtest_logpvalue, :mcc, :dips_tscore, :dips_sitescore, :correctclass_logp, :files
    def initialize(path)
      @path = path
      @files = {}
      File.open(path).each{|line|
        case(line)
        when /Discriminatory mutual information = (.*)/
          @mi = $1.to_f
        when /G-test = (.*)/
          @gtest_statistic = $1.to_f
        when /P\(Chi-Square\(G-Test\)\) = (.*)/
          @gtest_pvalue = $1.to_f
        when /Log-P\(Chi-Square\(G-Test\)\) = (.*)/
          @gtest_logpvalue = $1.to_f
        when /Matthews correlation coefficient = (.*)/
          @mcc = $1.to_f
        when /DIPS t-score = (.*)/
          @dips_tscore = $1.to_f
        when /DIPS site-score = (.*)/
          @dips_sitescore = $1.to_f
        when /log P correct classification = (.*)/
          @correctclass_logp= $1.to_f
        when /The log-likelihood for (.*) = (.*)/
          path = $1
          score = $2.to_f
          @files[path] = {} unless @files.key?(path)
          @files[$1][:log_likelihood] = score
        when /The expected posterior for (.*) = (.*)/
          path = $1
          score = $2.to_f
          @files[path] = {} unless @files.key?(path)
          @files[$1][:expected_posterior] = score
        when /The posterior for sequences with at least one motif of (.*) = (.*) \/ (.*) = (.*)/
          path = $1
          n = $2.to_f
          m = $3.to_f
          p = $4.to_f
          @files[path] = {} unless @files.key?(path)
          @files[path][:positive] = n
          @files[path][:total] = m
          @files[path][:p] = p
        end
      }
    end
    def to_s
      s = [@path, @mi, @gtest_statistic, @gtest_pvalue, @gtest_logpvalue, @mcc, @dips_tscore, @dips_sitescore, @correctclass_logp].map{|x| x.to_s }.join("\t")
      files.sort.each{|p, v|
        s += "\t#{v[:positive]}"
        s += "\t#{v[:total]}"
        s += "\t#{v[:p]}"
      }
      s
    end
  end
  def information_entropy
    e = 0
    (@n_states - @first_state).times{|i|
      state = @first_state + i
      @n_emissions.times{|j|
        e -= @emission[i][j] * Math::log(@emission[i][j]) if @emission[i][j] > 0
      }
    }
    e / Math::log(2.0)
  end
end

