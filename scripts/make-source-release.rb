#!/usr/bin/env ruby

require 'tempfile'

$noisy = false

if ARGV.length < 1
  puts "Please provide the output path."
  exit -1
end

out_path = ARGV[0]

unless out_path =~ /.tar.gz$/
  puts "Error: output path does not end in \".tar.gz\""
  exit -1
end

out_path = out_path

branch = "master"
`git branch`.each_line{|line|
  line.strip!
  if line =~ /\* (.*)/
    branch = $1
    if branch =~ /\(detached from (.*)\)/
      branch = $1
    end
  end
}

version = `git describe`.strip.gsub(/-/,".")
sha1 = `git rev-parse HEAD`.strip

if $noisy
  puts "branch = #{branch}"
  puts "version = #{version}"
  puts "sha1 = #{sha1}"
end

rel_path = "discrover-#{version}"

Dir.mktmpdir("discrover-source-archive"){|dir|
  tmp_archive_path = dir + "/archive.tar.gz"
  [
    "git archive #{branch} --prefix #{rel_path}/ | gzip > #{tmp_archive_path}",
    "cd #{dir} && tar -xzvf #{tmp_archive_path}",

    # hard-code some variables that would ordinarily be set during configuration using git routines
    "cd #{dir} && sed -i -e \"s/@GIT_DESCRIPTION@/#{version}/\" #{rel_path}/doc/discrover-manual.tex",
    "cd #{dir} && sed -i -e \"s/@GIT_DESCRIPTION@/#{version}/\" #{rel_path}/src/git_config.hpp.in",
    "cd #{dir} && sed -i -e \"s/@GIT_SHA1@/#{sha1}/\" #{rel_path}/src/git_config.hpp.in",
    "cd #{dir} && sed -i -e \"s/@GIT_BRANCH@/#{branch}/\" #{rel_path}/src/git_config.hpp.in",

    "cd #{dir} && tar -czvf #{out_path} #{rel_path}"
  ].each{|cmd|
    puts cmd if $noisy
    `#{cmd}`
  }
}
