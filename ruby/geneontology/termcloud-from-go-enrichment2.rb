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
    opts.on( '-f', 'use FDR < 0.05 instead of p < 0.05 filter'
      ){options[:filterfdr] = true}
    opts.on( '-c FILE', 'color mapping between filename and color code (without #)'
      ){|file| options[:colorfile] = file}
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
def map_color(ontology, direction)
  return '0065ab' if ontology == "BP" and direction == "O"
  return '93002d' if ontology == "MF" and direction == "O"
  return '1f9300' if ontology == "CC" and direction == "O"
  return '4d4d4d' if ontology == "BP" and direction == "U"
  return '4d4d4d' if ontology == "MF" and direction == "U"
  return '4d4d4d' if ontology == "CC" and direction == "U"
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
def load_color_map(file)
  str2color = Hash.new
  f = File.open(file, "r")
  while (line = f.gets)
    line.chomp!
    str, color = line.split("\t")[0,2]
    str2color[str] = color
  end
  return str2color
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
def parse_ORA_file(file, go2name, filterfdr=nil, color=nil)
  fw = File.open(file + ".termcloud", 'w')
  IO.foreach(file) do |line|
    direction, ontology, goid, p, fdr = line.split("\t")
    p, fdr = p.to_f, fdr.to_f
    next if p > 0.05
    next if filterfdr and fdr > 0.05
    p = '1e-200'.to_f if p == 0.0
    size = -1.0* Math.log(p)
    term = go2name[goid]
    STDERR.puts "could not find name for GO term #{goid}" unless term
    col = color ? color : map_color(ontology, direction)
    fw.puts [term, sprintf('%.2f', size), col].join(":")
  end
  fw.close
end

# =============================================================================
# === M A I N =================================================================
# =============================================================================

options = get_opt()
abort("directory does not exist - aborting.") unless File.exists?(options[:dir]) and File.directory?(options[:dir])
abort("go2name mapping file does not exist - aborting.") if options[:go2name] and not File.exists?(options[:go2name])
go2name = get_go2name(options[:go2name])
str2color = load_color_map(options[:colorfile]) if options[:colorfile]
Dir.glob(options[:dir] + '/*.ORA').each do |file|
  filename = File.basename(file)
  color = nil
  color = str2color[filename] if options[:colorfile]
  parse_ORA_file(file, go2name, options[:filterfdr], color)
end
