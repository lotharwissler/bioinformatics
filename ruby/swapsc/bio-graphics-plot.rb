#!/usr/bin/ruby
# generates a diagram of where in the sequence accelerated evolution / positive selection / negative selection took place
require 'rubygems'
require 'bio-graphics'

DEBUG = false



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

if ARGV[0] and not File.exists?(ARGV[0])
  puts "error: invalid path to file specified."
  ARGV[0] = nil
end

unless ARGV[0] or ARGV[1]
  puts "generates a diagram of where in the sequence accelerated evolution / positive selection / negative selection took place\n"
  puts "expected format [tab-delimited]:"
  puts "PANEL length"
  puts "TRACK name label"
  puts "FEATURE range color [label]"
  puts "usage: visualize-swapsc.rb <flatfile> <outfile>\n"
  exit 1
end

# === MAIN ====================================================================
# =============================================================================

# 1. read flatfile and save the input
# 2. process input, create the plot
panel = nil
track = nil
tracks = Array.new
features = Array.new

f = File.open( ARGV[0], "r" )
#STDERR.print( ARGV[0] + "\t" )
while line = f.gets
  next if line == nil
  line.chomp!
  cols = line.split("\t")
  if cols[0] == "PANEL"
    panel = Bio::Graphics::Panel.new( cols[1].to_i, :width => 800, :format => :png )
  elsif cols[0] == "TRACK"
    i, name, label = line.split("\t")
    if label == "true"
      label = true
    else
      label = false
    end
    track = panel.add_track(name, :label => label)
  elsif cols[0] == "FEATURE"
    if line.split("\t").length == 4
      i, range, color, label = line.split("\t")
      color = color.split(',').collect{|c| c.to_f}
      track.add_feature( Bio::Feature.new("feature", range), :colour => color, :label => label )
    else
      i, range, color = line.split("\t")
      color = color.split(',').collect{|c| c.to_f}
      track.add_feature( Bio::Feature.new("feature", range), :colour => color )
    end
  else
    STDERR.puts "unknown line descriptor \"#{cols[0]}\"" unless cols[0].nil?
  end
end
f.close
panel.draw(ARGV[1])

#STDERR.puts "done."
