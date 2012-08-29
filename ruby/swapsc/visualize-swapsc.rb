#!/usr/bin/ruby
# generates a diagram of where in the sequence accelerated evolution / positive selection / negative selection took place
require 'rubygems'
require 'bio-graphics'

DEBUG = false

$categories = { 
           "NS"  => { :color => [1,0,0], :legend => "Negative selection", :stats => 0 },
           "PS"  => { :color => [0.1,0.7,0.1], :legend => "Positive selection", :stats => 0 },
#           "HS"  => { :color => [0,0,0.7], :legend =>  "Hot sports", :stats => 0 },
#           "S"   => { :color => [0.5,0.5,0.5], :legend => "Saturation of synonymous sites", :stats => 0 },
#           "AdN+S" => { :color => [1,0.5,0], :legend => "Acceleration of non-synonymous substitutions + Saturation of synonymous sites", :stats => 0 },
           "AdN" => { :color => [0,0.7,1], :legend => "Acceleration of non-synonymous substitutions", :stats => 0 }
}

###############################################################################
class SwapscFeature

  attr_accessor :start, :stop, :category
  
  def initialize(start,stop,category)
    @start = start
    @stop = stop
    @category = category
    @added = false
  end

  def <=> other
    @start <=> other.start
  end  

  def added?
    return @added
  end

  def add_to_track(track)
    track.add_feature( Bio::Feature.new(@category, '%s..%s' % [ @start, @stop ]), :colour => $categories[@category][:color] )
    $categories[@category][:stats] += (@stop - @start +1)
    $categories[@category][:branchstats] += (@stop - @start +1)
    @added = true
  end
end
###############################################################################

###############################################################################
def correct_overlap( pf, cf )
  if cf.start <= pf.stop and cf.stop <= pf.stop # current included in previous
    if pf.category == 'NS'
      outerstop = +1 + pf.stop
      pf.stop = -1 + cf.start 
      nf = SwapscFeature.new(-1 + cf.stop,outerstop,"NS")
      return [pf, nf]
    elsif pf.category == 'PS' and cf.category == 'AdN'
      return [pf]
    elsif pf.category == 'AdN' and cf.category == 'PS'
      outerstop = +1 + pf.stop
      pf.stop = -1 + cf.start
      nf = SwapscFeature.new(-1 + cf.stop,outerstop,"AdN")
      return [pf, cf, nf]
    end
  else # current overlaps to the right (but is longer)
    if pf.category == 'NS' or cf.category == 'NS'
      newstop = -1 + cf.start
      newstart = +1 + pf.stop
      pf.stop = newstop
      cf.start = newstart
      return [pf, cf]
    elsif pf.category == 'AdN' and cf.category == 'PS'
      pf.stop = -1 + cf.start
      return [pf, cf]
    elsif pf.category == 'PS' and cf.category == 'AdN'
      cf.start = +1 + pf.stop
      return [pf, cf]
    end
  end
  return [pf, cf]
end
###############################################################################

if ARGV[0] and not File.exists?(ARGV[0])
  puts "error: invalid path to file specified."
  ARGV[0] = nil
end

unless ARGV[0]
  puts "generates a diagram of where in the sequence accelerated evolution / positive selection / negative selection took place\n"
  puts "usage: visualize-swapsc.rb <swapsc.out> [-r]"
  puts "\t-r\tremove overlapping features"
  puts "will create <swapsc.out>.png\n"
  exit 1
end

$annotation = true
tmpparts = ARGV[0].split('.')
annofilename = tmpparts[0..2].join('.') + '.ids.annotation'
#STDERR.puts "checking file #{annofilename} ..."
if File.exists?(annofilename)
  STDERR.puts "annotation file found: #{annofilename}"
else
  annofilename = '../' + annofilename
  #STDERR.puts "checking file #{annofilename} ..."
  if File.exists?(annofilename)
    STDERR.puts "annotation file found: #{annofilename}"
  else
    STDERR.puts "no annotation file found."
    $annotation = false
  end
end

$removeoverlaps = false
$removeoverlaps = true if ARGV[1] == '-r'

# === MAIN ====================================================================
# =============================================================================


globals = Hash.new
branches = Hash.new
branchsequences = Hash.new

f = File.open( ARGV[0], "r" )
  
  STDERR.print( ARGV[0] + "\t" )
  line = f.gets.chomp
  line.scan(/^Number of sequences =\s*(\d+)/) { |match| globals.store(:n, match.to_s.to_i) }
  #puts "Number of sequences: #{globals.fetch(:n)}"
  line = f.gets.chomp
  line.scan(/^Length of alignment \(nucleotides\) =\s*(\d+)/) { |match| globals.store(:length, match.to_s.to_i) }
  #puts "Length alignment: #{globals.fetch(:length)}"

  line = f.gets.chomp
  1.upto( globals.fetch(:n) ) do |i|
    line = f.gets.chomp
    branches.store(i, line)
    puts "branchpoint #{i}: #{line}" if DEBUG
    line = f.gets.chomp
    branchsequences.store(i, line)
    puts "  seq: #{line.slice(0,80)}..." if DEBUG
  end
  
  # Branches:
  # ----------
  line = f.gets.chomp while line !~ /^Branches:/
  line = f.gets.chomp
  line = f.gets.chomp
  while line =~ /^\d+\s+:\s+\d+\.{3}\d+/
    h = Hash.new
    line.scan(/^(\d+)/) { |match| h.store(:key, match.to_s.to_i) }
    line.scan(/^\d+\s+:\s+(\d+)\.{3}\d+/) { |match| h.store(:x, match.to_s.to_i) }
    line.scan(/^\d+\s+:\s+\d+\.{3}(\d+)/) { |match| h.store(:y, match.to_s.to_i) }
    #puts "#{h.fetch(:x)} #{h.fetch(:y)}"
    value = '(' + branches.fetch(h.fetch(:x)) + ',' + branches.fetch(h.fetch(:y)) + ')'
    puts "branch #{:key} #{h.fetch(:key)}: #{value} (#{h.fetch(:x)},#{h.fetch(:y)})" if DEBUG
    branches.store(h.fetch(:key), value)
    line = f.gets.chomp
  end

  # Ancestral sequences inferred by MP:
  # -----------------------------------
  line = f.gets.chomp while line !~ /^Ancestral sequences inferred by MP:/
  line = f.gets.chomp
  line = f.gets.chomp
  while line =~ /^node/
    descr, seq = line.split
    branch = (descr.match(/^node(\d+):/)[1]).to_s.to_i
    branchsequences.store(branch, seq)
    #puts "ancestral sequence #{branch} => #{seq.slice(0,80)}"
    line = f.gets.chomp
  end

  if DEBUG
    branches.each do |key,value| 
      puts "branch #{key} => #{value}"
      puts "  seq => #{branchsequences.fetch(key)}" if branchsequences.key?(key)
    end
  end

  line = f.gets.chomp while line !~ /^mean w =/
  line.scan(/^mean w =\s+(\S+);\s+/) { |match| globals.store(:omega, match.to_s) }
  #puts "omega: #{globals.fetch(:omega)}"
  
  
  panel = Bio::Graphics::Panel.new( globals.fetch(:length), :width => 800, :format => :png )
  # add annotation to diagram
  if $annotation
    annofile = File.new(annofilename, 'r')
    annotationtext = annofile.readline.chomp
    annotationtext = annotationtext.gsub('"','')
    annofile.close
    STDERR.puts "annotation line: #{annotationtext}"
    track = panel.add_track('TAIR annotation', :label => true, :colour => [1.0,0.8,0.3])
    track.add_feature( Bio::Feature.new('tair8', '%s..%s' % [ 0, globals.fetch(:length) ]), :label => annotationtext )
  end

  # GET LIST OF SIGNALS
  line = f.gets.chomp while line !~ /^={20,}/

  while line !~ /^Proportion of codon sites under selective constraints/ # read all branches
    while line !~ /^Proportion of codon sites under selective constraints/ and line !~ /\d+\.{2}\d+/ # read single branch
      # puts "skipping line: #{line}"
      line = f.gets.strip # while line !~ /\d+\.{2}\d+/
    end
    if line =~ /\d+\.{2}\d+/
      puts "positions: #{line}" if DEBUG
      bs = line.split('..')
      x = branches.fetch(bs[0].to_i)
      y = branches.fetch(bs[1].to_i)
      name = "#{x} : #{y}"
      track = panel.add_track(name, :label => false, :colour => [0,0,0])
      # display gaps
      if branchsequences.key?(bs[0].to_i) and branchsequences.key?(bs[1].to_i)
        i = 0
        seqx = branchsequences.fetch(bs[0].to_i)
        seqy = branchsequences.fetch(bs[1].to_i)
        puts "seqx (#{bs[0]}): #{seqx}" if DEBUG
        puts "seqy (#{bs[1]}): #{seqy}" if DEBUG
        start, stop = nil, nil
        while i < seqx.length do
          if (seqx[i,1] == '-' or seqy[i,1] == '-') 
            if start.nil?
              start = i 
            end
            stop = i

          else
            if not start.nil? and not stop.nil?
               track.add_feature( Bio::Feature.new("gap", '%s..%s' % [ start, stop ]), :colour => [0.55,0.55,0.55] )
               puts "added gap in #{name} between #{start} and #{stop}" if DEBUG
               start, stop = nil, nil
            end
          end
          i += 1
        end
      else
        STDERR.puts "no sequence pair found for #{name}"
      end
      # /display gaps
      line = f.gets.chomp
      features = Hash.new
      while line.split.size >= 10
        columns = line.split
        significance = columns.slice(6,3).to_s
        unless significance == "P>0.05"
          positions = columns[0]
          category = columns.slice(9,columns.size).to_s
          if $categories.has_key?(category)
            puts "add feature #{positions} \"#{category}\"" if DEBUG
            features[category] = Hash.new unless features.has_key?(category)
            start, stop = positions.split('..')
            features[category][start.to_i] = stop.to_i 
          end
        end
        line = f.gets.chomp
      end
      
      sortedfeatures = Array.new
      features.each do |category, positionhash| # iterate features
        next if positionhash.size == 0
        sortedpositions = positionhash.sort
        fsf = nil
        for e in sortedpositions # iterate positions
          nstart, nstop = e
          puts "#{category} #{nstart} #{nstop}" if DEBUG
          unless fsf
            fsf = SwapscFeature.new(nstart,nstop,category)
            next
          end
          if nstart <= fsf.stop
            fsf.stop = nstop
          else
            sortedfeatures << fsf ? $removeoverlaps : fsf.add_to_track(track) 
            fsf = SwapscFeature.new(nstart,nstop,category)
          end
        end # /iterate positions
        
        sortedfeatures << fsf ? $removeoverlaps : fsf.add_to_track(track) unless fsf.added?

      end # /iterate features

      if $removeoverlaps 
          sortedfeatures.sort!
          sortedfeatures.each_index do |index|
              unless index == 0
                currentfeature = sortedfeatures[index]
                prevfeature = sortedfeatures[(index -1)]
                if currentfeature.start < prevfeature.stop:
                  #puts sortedfeatures.inspect
                  #STDOUT.puts "\n#{name}\toverlap found:\n"
                  #STDOUT.puts "\t%s\tstart: %s\tstop:\t%s" % [prevfeature.category, prevfeature.start, prevfeature.stop]
                  #STDOUT.puts "\t%s\tstart: %s\tstop:\t%s" % [currentfeature.category, currentfeature.start, currentfeature.stop]
                  resarray = correct_overlap(prevfeature,currentfeature)
                  sortedfeatures[index -1] = nil
                  sortedfeatures[index] = nil
                  resarray.each {|r| sortedfeatures << r }
                  sortedfeatures.compact!
                  sortedfeatures.sort!
                  retry
                end
              end
          end # /sortedfeatures.each
          sortedfeatures.compact!
          $categories["PS"][:branchstats] = 0
          $categories["NS"][:branchstats] = 0
          $categories["AdN"][:branchstats] = 0
          sortedfeatures.each do |feat| 
            feat.add_to_track(track) 
          end
          printlist = Array.new
          printlist << ARGV[0]
          printlist << name
          $categories.each { |key,hash| printlist << $categories[key][:branchstats] }
          STDOUT.puts( printlist.join("\t") )
      end

    end # /read single branch 


  end # /read all branches


#  exit 3


  # LEGEND + STATS
  track = panel.add_track("Legend", :label => true)
  startpos = (globals.fetch(:length).to_f * 0.02).to_i
  endpos = (globals.fetch(:length).to_f * 0.98).to_i
  negative, positive = 0, 0
  #STDOUT.print "#{ARGV[0]}"
  $categories.each do |abbrv,hash|
    color = hash[:color]
    legend = hash[:legend]
    puts "legend: #{abbrv} #{legend} #{color}" if DEBUG
    stats = hash[:stats].to_f * 100 / (branches.size - globals.fetch(:n)) / globals.fetch(:length)
    if ["S","NS"].include?(abbrv)  
      negative += stats
    else
      positive += stats
    end
    stats = format("%.2f",stats)
    #STDOUT.print "\t#{abbrv}\t#{stats}"
    track.add_feature( Bio::Feature.new(abbrv, '%s..%s' % [ startpos, endpos ]), :colour => color, :label => "#{abbrv} (#{legend}): #{stats} %" )
  end
  positive = format("%.2f",positive)
  negative = format("%.2f",negative)
  #STDOUT.print "\t#{positive}\t#{negative}\n"

  panel.draw(ARGV[0] + ".gaps.png")
  STDERR.print( "done.\n" )

f.close
