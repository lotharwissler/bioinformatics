#!/usr/bin/ruby -w
# == Synopsis
# creates input for wordle.net/advanced to create a term cloud based on the six
# output files produced by the go-enrichment.py script
#

require 'optparse'
require 'rubygems'
require 'faster_csv'


FILES = ['topGO.over.Sig.CC.csv', 'topGO.over.Sig.BP.csv', 'topGO.over.Sig.MF.csv', 'topGO.under.Sig.CC.csv', 'topGO.under.Sig.BP.csv', 'topGO.under.Sig.MF.csv']

# =============================================================================
def get_opt
  options = Hash.new
  optparse = OptionParser.new do |opts|
    opts.banner = "Usage: #{$0} -d <dir>"
    options[:dir] = nil
    options[:go2name] = nil
    opts.on( '-d DIR', 'directory that contains the topGO.over.Sig*.csv files'
      ){|dir| options[:dir] = dir}
    opts.on( '-g FILE', 'gene ontology id to name mapping file, tab delimited, to look up shortened term names'
      ){|file| options[:go2name] = file}
  end
  begin
    optparse.parse!
    mandatory = [:dir]
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
def map_color(ontology, direction)
  return '0065ab' if ontology == "BP" and direction == "over"
  return '93002d' if ontology == "MF" and direction == "over"
  return '1f9300' if ontology == "CC" and direction == "over"
  return '4d4d4d' if ontology == "BP" and direction == "under"
  return '4d4d4d' if ontology == "MF" and direction == "under"
  return '4d4d4d' if ontology == "CC" and direction == "under"
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
def parse_topgo_file(file, go2name)
  direction, ontology = file.split('.')[1], file.split('.')[3]
  csvrows = FasterCSV.read(file, :quote_char => '"', :col_sep => ',', :row_sep => :auto, :headers => true)
  csvrows.each do |row|
    pfilter = row[10].to_f
    next unless pfilter < 0.05
    size = -1.0* Math.log(pfilter)
    term = row[2]
    term = go2name[row[1]] if go2name and term[-3,3] == '...'
    color = map_color(ontology, direction)
    puts [term, size.to_s, color].join(":")
  end
end

# =============================================================================
# === M A I N =================================================================
# =============================================================================

options = get_opt()
abort("directory does not exist - aborting.") unless File.exists?(options[:dir]) and File.directory?(options[:dir])
abort("go2name mapping file does not exist - aborting.") if options[:go2name] and not File.exists?(options[:go2name])
go2name = nil
go2name = get_go2name(options[:go2name]) if options[:go2name]
Dir.chdir(options[:dir])
FILES.each do |file|
  unless File.exists?(file)
    STDERR.puts "could not find file #{file} in the given dir. skipping..." unless File.exists?(file)
    next
  end
  parse_topgo_file(file, go2name)
end
