#---------------------------------------------------
require 'uri'
require 'bio'
require 'open-uri'
require 'net/http'

#Initialization hash - comment or remove parameters which are not going to be modified
params=
{
  :ALIGNMENTS =>'10',
  #:ALIGNMENT_VIEW =>'',
  #:AUTO_FORMAT =>'',
  #:COMPOSITION_BASED_STATISTICS =>'',
  :DATABASE =>'nr',
  #:DB_GENETIC_CODE =>'',
  :DESCRIPTIONS =>'10',
  #:ENDPOINTS =>'',
  #:ENTREZ_LINKS_NEW_WINDOW =>'',
  #:ENTREZ_QUERY =>'',
  #:EXPECT =>'',
  #:EXPECT_LOW =>'',
  #:EXPECT_HIGH =>'',
  :FILTER =>'L',
  #:FORMAT_ENTREZ_QUERY =>'',
  #:FORMAT_OBJECT =>'',
  #:FORMAT_TYPE =>'',
  #:GAPCOSTS =>'',
  #:GENETIC_CODE =>'',
  #:HITLIST_SIZE =>'',
  #:I_THRESH =>'',
  #:LAYOUT =>'',
  #:LCASE_MASK =>'',
  #:MEGABLAST =>'',
  #:MATRIX_NAME =>'',
  #:NCBI_GI =>'',
  #:NUCL_PENALTY =>'',
  #:NUCL_REWARD =>'',
  #:OTHER_ADVANCED =>'',
  #:PAGE =>'',
  #:PERC_IDENT =>'',
  #:PHI_PATTERN =>'',
  #:PROGRAM =>nil,
  #:PSSM =>'',
  #:QUERY_BELIEVE_DEFLINE =>'yes',
  #:QUERY_FROM =>'',
  #:QUERY_TO =>'',
  #:RID =>'',
  #:RESULTS_FILE =>'',
  #:SEARCHSP_EFF =>'',
  #:SERVICE =>'',
  #:SHOW_OVERVIEW =>'',
  #:THRESHOLD =>'',
  #:UNGAPPED_ALIGNMENT =>'',
  #:WORD_SIZE =>'',
  #:TARGET =>'',
}

#Example sequence, maybe in FASTA format or in plain text
#White spaces and not allowed signs will be removed
seq1='------------------'
sequence='>kjljljjl
ACDCrkrkDCDCD---xCDCDCDCD'
seq = ">gi|gggghghg|embl|mnm,mnbjbkbkbjkbkkjb
ATGcagctctttgtccgcgcccaggagctacacaccttcgaggtgaccggccaggaaacg
gtcgcccagatcaaggctcatgtagcctcactggagggcattgccccggaagatcaagtc
gtgctcctggcaggcgcgcccctggaggatgaggccactctgggccagtgcggggtggag
gccctgactaccctggaagtagcaggccgcatgcttggaggtaaagtccatggttccctg
gcccgtgctggaaaagtgagaggtcagactcctaaggtggccaaacaggagaagaagaag
aagaagacaggtcgggctaagcggcggatgcagtacaaccggcgctttgtcaacgttgtg
cccacctttggcaagaagaagggccccaatgccaactcttaa"

#-----------------------------------------------------------------------------------
# Versatile remote NCBI BLAST - built over BioRuby
# Each instance of the class represents single BLAST request
# Initialization defines the parameters of the request provided in the form of a hash

class NCBI_BLAST
  include Bio
  #BEGIN--------CONSTANTS-------------------------------------------------
  PROGRAMS=['blastn','tblastn','tblastx','blastp','blastx','']
  REGEX_AA=Regexp.new('[A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|U|V|W|X]*')
  REGEX_NA=Regexp.new('[a|c|g|t|u|>|\n]*')
  URL=URI.parse('http://www.ncbi.nlm.nih.gov/blast/Blast.cgi')
  #END----------CONSTANTS-------------------------------------------------
  
  #BEGIN--------ATTRIBUTES-------------------------------------------------
  attr :error_message
  attr :seq_type
  attr :header
  attr :rid
  attr :done
  @program
  @query
  @error_message=''
   
  #BEGIN--------INITIALIZE-------------------------------------------------
  def initialize(sequence, params={}, seq_type_detect=false)
    @query_params=
    {
      :PROGRAM => 'blastn'
    }.merge(params)
        
    error=false
    @done=false
    # Reads QUERY sequence and checks its validity. If sequence (in FASTA or plain format) contains
    # not allowed symbols removes them and autodetects the sequence type.
    # The header of the sequence in FASTA format is stored in @header variable (otherwise @header is nil)
    @header=sequence.slice(/^>(.*)/).to_s

    #seq_temp=sequence.delete sequence.slice(/^>(.*)/).to_s

    @query=Bio::Sequence.auto(sequence)
    # puts @query.seq.class.to_s
    case @query.seq.class.to_s
      when "Bio::Sequence::AA" then
      #@query=sequence.upcase.scan(REGEX_AA).to_s
      @seq_type='AA'
      when "Bio::Sequence::NA" then
       #@query=sequence.downcase.scan(REGEX_NA).to_s
       @seq_type='NA'
    end
    @query=sequence.to_s
    # Assigns BLAST PROGRAM -
    # If no PROGRAM provided 'blastn' is being tried to be assigned
    # If PROGRAM does not apply to the QUERY type the error is returned
    @program=@query_params[:PROGRAM]
    if not PROGRAMS.include?(@program)
      @error_message='UNKNOWN PROGRAM'
      error=true
    elsif ((@seq_type=='NA' and not PROGRAMS[0..2].include?(@program)) or (@seq_type=='AA' and not PROGRAMS[3..4].include?(@program)))
      @error_message= "\'%s\' SEQUENCE PROVIDED FOR \'%s\' PROGRAM!" % [@seq_type, @program.upcase]
      error=true
    end
    
    # Reports error /for debugging only?/
    if error
      puts "Problem with initialization: " + @error_message
    end
  end
  #END----------INITIALIZE-------------------------------------------------

  #BEGIN--------DO_BLAST-------------------------------------------------
  # Returns String value representing assigned RID
  # In case of the error returns nil
  
  def do_blast
    sub=@query_params.merge({'CMD'=>'Put', 'QUERY'=>@query})
    http = Net::HTTP::post_form(URL, sub)
    @rid=http.body.scan(/RID = (.*?)$/).to_s
  end
  #END----------DO_BLAST------------------------------------------------------

  #BEGIN--------FETCH_RESULTS-------------------------------------------------
  # Fetches the results of the commited blast search
  # Returns result of the search in XML format if successfully retrieved the data or false if error occurred
  # Accepts one optional argument:
  #                 true  - waits until BLAST search is finished and returns
  #                 false - checks if BLAST 
  def fetch_results(wait=true)
    if @rid!=nil  
      error=false
      res='http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=' + @rid
      while open(res).read.scan(/Status=(.*?)$/).to_s=='WAITING'
        wait == false ? break : puts("Status=WAITING")
        sleep(1)
      end  
      case open(res).read.scan(/Status=(.*?)$/).to_s
        when 'WAITING' then 
          error=true
          @error_message="STILL WAITING FOR RESULTS"
        when 'FAILED' then 
          error=true
          @error_message="Fetching of the data failed"
        when 'UNKNOWN' then
          error=true
          @error_message="RID UNKNOWN OR EXPIRED"
      end
      if not error
        res='http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID=' + rid.to_s
        @error_message =''
        @done=true
        @answer=open(res).read
      else return false
      end
    else 
      @error_message = "BLAST SEARCH DID NOT INITIATED."
      return false 
    end
  end
  #END----------FETCH_RESULTS------------------------------------------------
end
#------END of NCBI_BLAST CLASS DEFINITION------------------------------------

h=NCBI_BLAST.new(seq, params)
puts h.seq_type
h.do_blast
while h.fetch_results==false and h.error_message=="STILL WAITING FOR RESULTS"
  puts h.error_message
end
if h.done
#  puts h.fetch_results(false)
rep=Bio::Blast::Report.new(h.fetch_results)
  puts "Hits for " + rep.query_def + " against " + rep.db
  rep.each do |hit|
    hit.each do |hsp|
    print hit.target_id, "\t", hit.evalue, "\t", hit.definition,"\n"
    print "\t\t\t", hsp.evalue, "\t", hsp.identity, "\t", hsp.percent_identity(), "\n"
  end
end
else puts "ERROR"
end  
puts h.error_message
#h.process_output