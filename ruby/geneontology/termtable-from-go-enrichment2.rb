#!/usr/bin/ruby -w
# == Synopsis
# creates input for wordle.net/advanced to create a term cloud based on the six
# output files produced by the go-enrichment.py script
#

require 'optparse'
require 'rubygems'
require 'faster_csv'
#require 'rsruby'


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
    opts.on( '-m MIN', 'minimum number of species in which a term has to be found significant'
      ){|s| options[:min] = s.to_i}
    opts.on( '-f', 'use FDR < 0.05 instead of p < 0.05 filter'
      ){options[:filterfdr] = true}
    opts.on( '-o', 'evaluate over-representation' 
      ){options[:over] = true} 
#    opts.on( '-u', 'evaluate under-representation' 
#      ){options[:under] = true} 
  end
  begin
    optparse.parse!
    mandatory = [:dir, :go2name, :min]
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
def parse_ORA_file(file, filterfdr, doOver, doUnder)
  terms = Hash.new
  IO.foreach(file) do |line|
    direction, ontology, goid, p, fdr = line.split("\t")
    next if direction == 'O' and not doOver
    next if direction == 'U' and not doUnder
    p, fdr = p.to_f, fdr.to_f
    next if p > 0.05
    next if filterfdr and fdr > 0.05
    terms[goid] = filterfdr ? fdr : p
  end
  return terms
end

# =============================================================================
# === M A I N =================================================================
# =============================================================================

options = get_opt()
abort("directory does not exist - aborting.") unless File.exists?(options[:dir]) and File.directory?(options[:dir])
abort("go2name mapping file does not exist - aborting.") if options[:go2name] and not File.exists?(options[:go2name])
go2name = get_go2name(options[:go2name])
enrichedHash = Hash.new
Dir.glob(options[:dir] + '/*.ORA').each do |file|
  key = File.basename(file).split('.').first
  enrichedHash[key] = parse_ORA_file(file, options[:filterfdr], options[:over], options[:under])
end

species = enrichedHash.keys
STDOUT.puts((["GO.ID", "GO.NAME"] + species.collect{|s| s.upcase}).join("\t"))
allterms = enrichedHash.values.collect{|v| v.keys}.inject{|union, array| union + array}.uniq
allterms.each do |goid|
  occurrence = enrichedHash.select{|s, gohash| gohash.key? goid}.count
  next if occurrence < options[:min]
  STDOUT.print [goid, go2name[goid]].join("\t")
  out = species.collect{|s| (enrichedHash[s].key?(goid) ? Math.log(-1*Math.log(enrichedHash[s][goid])).to_s : "0")}
  STDOUT.puts "\t" + out.join("\t")
end
STDOUT.close

#r = RSRuby.instance
#r.assign('d', r.read.csv(options[:dir] + '/termtable.tab', :header => true, :sep => "\t"))
#r.assign('m', as.matrix(r.d[,3:species.count+2]))
#r.rownames(r.d) = r.d[,2]
#r.library('gplots')
#r.library('RColorBrewer')
#r.pdf(options[:dir] + '/termtable.pdf')
#r.heatmap.2(r.m, :col => brewer.pal(3, "Blues"))
#r.eval_R("dev.off()")
