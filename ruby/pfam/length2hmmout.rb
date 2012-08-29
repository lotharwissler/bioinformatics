#!/usr/bin/ruby -w
HEADER = /^>(\S+)\s?.?/
unless (ARGV.size == 2)
	puts "Usage: #{$0} fasta hmmout [NOTE: will change input hmmout!]"
	exit
end
lengths	= Hash.new
seq 		= String.new
pid			= nil 
f 			= File.open(ARGV[0], "r")
c				= 0
while(line = f.gets)
	line.chomp!
	if (m = HEADER.match(line))
		lengths[pid] = seq.length.to_s unless (pid.nil?)
		pid = m[1]
		seq = String.new
		c += 1
		STDERR.print "\r*** Reading fasta entries: #{c}... "
		next
	end
	seq += line	
end
lengths[pid] = seq.length unless (pid.nil?)
f.close
STDERR.puts "done."
oldhmmout = Array.new
IO.foreach(ARGV[1]) {|x| oldhmmout << x}
f = File.open(ARGV[1], "w")
oldhmmout.each do |line|
	next if (/^#.+/.match(line))
	fields = line.split
	unless (lengths.has_key?(fields[0]))
			puts "*** NO LENGTH FOUND FOR >#{fields[0]}<"
			present = false
			next
	end
	line.chomp!
	f.puts lengths[fields[0]].to_s + "\t" + line + "\n"
end
f.close
