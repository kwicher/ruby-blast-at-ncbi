#---------------------------------------------------
require 'uri'
require 'bio'
require 'open-uri'
require 'net/http'

#Initialization hash
params=
{
  :PROGRAM =>'',
  :ALIGNMENTS =>'',
  :ALIGNMENTS =>'',
  :ALIGNMENT_VIEW =>'',
  :AUTO_FORMAT =>'',
  :CMD =>'',
  :COMPOSITION_BASED_STATISTICS =>'',
  :DATABASE =>'',
  :DB_GENETIC_CODE =>'',
  :DESCRIPTIONS =>'',
  :ENDPOINTS =>'',
  :ENTREZ_LINKS_NEW_WINDOW =>'',
  :ENTREZ_QUERY =>'',
  :EXPECT =>'',
  :EXPECT_LOW =>'',
  :EXPECT_HIGH =>'',
  :FILTER =>'',
  :FORMAT_ENTREZ_QUERY =>'',
  :FORMAT_OBJECT =>'',
  :FORMAT_TYPE =>'',
  :GAPCOSTS =>'',
  :GENETIC_CODE =>'',
  :HITLIST_SIZE =>'',
  :I_THRESH =>'',
  :LAYOUT =>'',
  :LCASE_MASK =>'',
  :MEGABLAST =>'',
  :MATRIX_NAME =>'',
  :NCBI_GI =>'',
  :NUCL_PENALTY =>'',
  :NUCL_REWARD =>'',
  :OTHER_ADVANCED =>'',
  :PAGE =>'',
  :PERC_IDENT =>'',
  :PHI_PATTERN =>'',
  :PROGRAM =>'',
  :PSSM =>'',
  :QUERY_BELIEVE_DEFLINE =>'',
  :QUERY_FROM =>'',
  :QUERY_TO =>'',
  :RID =>'',
  :RESULTS_FILE =>'',
  :SEARCHSP_EFF =>'',
  :SERVICE =>'',
  :SHOW_OVERVIEW =>'',
  :THRESHOLD =>'',
  :UNGAPPED_ALIGNMENT =>'',
  :WORD_SIZE =>'',
  :TARGET =>'',
}

#Example sequence
sequence='ACDCDCDCDCD---xCDCDCDCD'

#Versatile remote NCBI BLAST - built over BioRuby
#Each instance of the class represents single BLAST request
#Initialization defines the parameters of the request provided in the form of a hash
class NCBI_BLAST
  include Bio
  PROGRAMS=['blastn','blastp','blastx','tblastn','tblastx','']
  attr :program
  attr :query
  def initialize(sequence, params={})
    query_params=
    {
      :PROGRAM => 'blastn'
    }.merge(params)
   
   @success=true
   
   #Reads QUERY sequence and checks its validity. If QUERY contains not allowed symbols
   #removes them and autodetects its type
   seq_temp=sequence
   @query=Bio::Sequence.auto(sequence)
   case @query.seq.class.to_s
     when "Bio::Sequence::AA": puts seq_temp1=seq_temp.scan(/[A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|U|V|W|X]*/).to_s + " AA"
     when "Bio::Sequence::NA": puts seq_temp1=seq_temp.scan(/[A|C|G|T|U]*/).to_s + " NA"
    end 
    
   #Assign BLAST PROGRAM - if empty, 'blastn' is used
   @program = params[:PROGRAM]=='' ? 'blastn' : params[:PROGRAM]
   #puts @program
    if not PROGRAMS.include?(@program)
      @success='PROGRAM'
    end
    if(@success != true) 
      puts "Problem with initialization - check \"" + @success +"\" parameter"
    end
  end
end
h=NCBI_BLAST.new(sequence, params)
puts h.class
puts h.program
puts h.query
seq = Bio::Sequence::NA.new("
atgcagctctttgtccgcgcccaggagctacacaccttcgaggtgaccggccaggaaacg
gtcgcccagatcaaggctcatgtagcctcactggagggcattgccccggaagatcaagtc
gtgctcctggcaggcgcgcccctggaggatgaggccactctgggccagtgcggggtggag
gccctgactaccctggaagtagcaggccgcatgcttggaggtaaagtccatggttccctg
gcccgtgctggaaaagtgagaggtcagactcctaaggtggccaaacaggagaagaagaag
aagaagacaggtcgggctaagcggcggatgcagtacaaccggcgctttgtcaacgttgtg
cccacctttggcaagaagaagggccccaatgccaactcttaa")
data={'CMD'=>'Put', 'QUERY'=>seq, 'NCBI_GI'=>'yes', 'PROGRAM'=>'blastn', 'FILTER'=>'L', 'DATABASE'=>'nr', 'ALIGNMENTS'=>'100', 'DESCRIPTIONS'=>'100'}
url=URI.parse('http://www.ncbi.nlm.nih.gov/blast/Blast.cgi')
#http = Net::HTTP::post_form(url, data)
#rid=http.body.scan(/RID = (.*?)$/)
#res='http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=' + rid.to_s
puts seq.to_s.scan(/[^>]$.(.*)/m)
#puts res
=begin answer=open(res).read
while open(res).read.scan(/Status=(.*?)$/).to_s=='WAITING' do
    sleep(1)
    puts "Status=WAITING"
   end  
puts "Status=" + open(res).read.scan(/Status=(.*?)$/).to_s
res='http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID=' + rid.to_s
answer=open(res).read 
puts "OK:\n"
hits=answer.scan(/<Hit_id>(.*?)<\/Hit_id>/)
hits.each{|hit| puts hit}
=end
