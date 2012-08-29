#/usr/bin/ruby
=begin
=end

class GOterm
  attr_accessor :id, :name, :namespace, :parents
  def initialize
    @parents = Array.new
  end
end

def load_obo_definition(file)
  goterm = Hash.new
  obofile = File.open(file)
  while line = obofile.gets.chomp
    if line =~ /^\[Term\]/
      g = GOterm.new
    elsif line =~ /^id:/
      g.id = line.scan(/^id:\s+(GO:\d+)/).first.first
    elsif line =~ /^name:/
      g.name = line.scan(/^name:\s+(.*)$/).first.first
    elsif line =~ /^namespace:/
      g.namespace = line.scan(/^namespace:\s+(\S+)$/).first.first
    elsif line =~ /^is_a:/
      g.parents << line.scan(/^is_a:\s+(GO:\d+)/).first.first
    elsif line =~ /^$/
      goterm[g.id] = g
    end
  end
  return goterm
end



