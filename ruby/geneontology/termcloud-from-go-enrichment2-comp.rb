#!/usr/bin/ruby -w
# == Synopsis
# creates input for wordle.net/advanced to create a term cloud based on the six
# output files produced by the go-enrichment.py script
#

require 'optparse'
require 'rubygems'
require 'faster_csv'


# =============================================================================
def get_opt
  options = Hash.new
  optparse = OptionParser.new do |opts|
    opts.banner = "Usage: #{$0} -d <dir>"
    options[:dir] = nil
    options[:go2name] = nil
    opts.on( '-d DIR', 'directory that contains *.ORA files produced with go-enrichment2.py'
      ){|dir| options[:dir] = dir}
    opts.on( '-g FILE', 'gene ontology id to name mapping file, tab delimited, to look up shortened term names'
      ){|file| options[:go2name] = file}
  end
  begin
    optparse.parse!
    mandatory = [:dir, :go2name]
    missing = mandatory.select{|param| options[param].nil?}
    if not missing.empty?
      puts "Missing options: #{missing.join(', ')}"
      puts optparse 
      exit
    end
  rescue OptionParser::InvalidOption, OptionParser::MissingArgument
    puts $!.to_s
    puts optparse
    exit
  end
  return options
end


# modify to adjust colors
# =============================================================================
def map_color(species)
  return '004d84' if species == "Acep"
  return '5abbff' if species == "Aech"
  return 'b4b4b4' if species == "Amel"
  return '000000' if species == "Hsal"
  return '93002d' if species == "Cflo"
  return 'ff6d9a' if species == "Lhum"
  return '145e00' if species == "Pbar"
  return '71f84c' if species == "Sinv"
  return '803300' if species == "Nvit"
  return 'ff9955' if species == "Dmel"
end

# =============================================================================
def statusbar(progress, message="", width=40)
  progressbar = "=" * (progress*width).to_i
  progressbar << " " while progressbar.length < width
  STDERR.print "\r   0% #{progressbar} 100% "
  STDERR.print "[#{message}]" unless message.empty?
  STDERR.print "\n" if progress == 1.0 
end

# =============================================================================
def get_go2name(file)
  go2name = Hash.new
  f = File.open(file, "r")
  while (line = f.gets)
    line.chomp!
    id, name = line.split("\t")[0,2]
    go2name[id] = name
  end
  return go2name
end

# =============================================================================
def parse_ORA_file(file, go2name)
  species = file.split(".")[-2]
  IO.foreach(file) do |line|
    direction, ontology, goid, p, fdr = line.split("\t")
    p, fdr = p.to_f, fdr.to_f
    next if fdr > 0.05
    size = -1.0* Math.log(p)
    term = go2name[goid]
    STDERR.puts "no name found for #{goid}" unless term
    color = map_color(species)
    STDOUT.puts [term, sprintf('%.2f', size), color].join(":")
  end
end

# =============================================================================
# === M A I N =================================================================
# =============================================================================

options = get_opt()
abort("directory does not exist - aborting.") unless File.exists?(options[:dir]) and File.directory?(options[:dir])
abort("go2name mapping file does not exist - aborting.") if options[:go2name] and not File.exists?(options[:go2name])
go2name = get_go2name(options[:go2name])
Dir.glob(options[:dir] + '/*.ORA').each do |file|
  parse_ORA_file(file, go2name)
end
