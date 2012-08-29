#!/usr/bin/ruby

require 'optparse'

PFAMFILE = "/global/databases/pfam/current/pfam_scan_db/Pfam-A.hmm"

class String
  def valid_float?
    # The double negation turns this into an actual boolean true - if you're 
    # okay with "truthy" values (like 0.0), you can remove it.
    !!Float(self) rescue false
  end
end


# =============================================================================
def get_opt
  options = {}
  optparse = OptionParser.new do |opts|
    opts.banner = "Usage: #{$0} -f <file> -c <value>"
    opts.on( '-f FILE or DIR', 'single hmmout file (pfam_scan output with first column = protein length), or a directory where all *.hmmout files will be processed' 
      ){|file| options[:hmmfile] = file}
    opts.on( '-c CUTOFF', '[evalueFloat|GA|TC|NC]'
      ){|v| options[:cutoff] = v}
  end 
  begin
    optparse.parse!
    mandatory = [:hmmfile, :cutoff]
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

def get_cutoffs(file=PFAMFILE)
  cutoffHash = Hash.new
  capture = %w( NAME GA NC TC )
  @name = nil
  reader = File.open(file, 'r')
  while (line = reader.gets)
    entry = {} if line[0,6] == 'HMMER3'
    capture.each{|e| entry[e] = line.split[1] if line[0,e.length] == e }
    if line[0,2] == "//" 
      if entry.length != capture.count
        STDERR.puts "FATAL ERROR: not all required fields found for an entry: #{entry.inspect}"
        next
      end
      cutoffHash[entry['NAME']] = entry
    end
  end
  return cutoffHash
end


# ==============================================================================
def filter_hmmout(file, cutoff)
  fw = File.open(file + "." + cutoff, 'w')
  f = File.open(file, 'r')
  if cutoff.valid_float? # e-value cutoff given
    e = cutoff.to_f
    f.each{|line| cols = line.chomp.split; fw.puts cols.join("\t") if cols[13].to_f < e}
  else 
    e = cutoff if ['GA', 'TC', 'NC'].include?(cutoff)
    abort("invalid value given for cutoff method (#{cutoff}). allowed values are GA, NC, and TC.") if e.nil?
    cutoffHash = get_cutoffs()
    puts "--- cutoffHash: #{cutoffHash.count} ---"
    f.each{|line| 
      cols = line.chomp.split;
      name, bitscore = cols[7], cols[12].to_f
      puts name, cutoffHash[name]
      next unless name[0,6] == 'Pfam-B' or bitscore > cutoffHash[name][e].to_f
      fw.puts cols.join("\t") 
    }
  end
  f.close
  fw.close
end


# ==============================================================================
# =MAIN=========================================================================
# ==============================================================================

options = get_opt()
unless File.directory?(options[:hmmfile])
  filter_hmmout(options[:hmmfile], options[:cutoff])
else
  Dir.glob(options[:hmmfile] + '/*.hmmout').each{|hmmfile| filter_hmmout(hmmfile, options[:cutoff])}
end

