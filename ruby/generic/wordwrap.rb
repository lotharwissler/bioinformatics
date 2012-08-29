#!/usr/bin/ruby
require 'ftools'

def exit_with_usage
  STDOUT.print """ 
  #{$0} <file> [<line length>]
  
  this script inserts newlines in front of all words that would exceed
  a given threshold for max line length.
  default max length: 80

  """
  exit(1)
end

exit_with_usage unless ARGV.length > 0
exit_with_usage unless File.exists? ARGV[0]
MAXLENGTH = (ARGV[1] || 80).to_i

STDERR.puts "INPUT FILE:\t%s" % ARGV[0]
STDERR.puts "MAXLENGTH:\t%s" % MAXLENGTH

f = File.open(ARGV[0])
while line = f.gets
  if line.length < 80
    STDOUT.print line
  else
    words = line.chomp.split
    first = line[0..0]
    pos = 0
    newline = Array.new
    words.each do |word|
      newline << word
      if newline.join(" ").length > MAXLENGTH
        newline[-1] = "\n"
        STDOUT.print newline.join(" ")
        if first == "#" or first == "%"
          newline = [first, word]
        else
          newline = [word]
        end
      end
    end
    STDOUT.puts newline.join(" ").chomp
  end
end
f.close

STDERR.puts "STATUS:   \tdone.\n"
