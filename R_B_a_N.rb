#---------------------------------------------------
require 'uri'
require 'bio'
require 'open-uri'
require 'net/http'

#Initialization hash - comment or remove parameters which are not going to be modified
params=
{
  :ALIGNMENTS =>'1',
  #:ALIGNMENT_VIEW =>'',
  #:AUTO_FORMAT =>'',
  #:COMPOSITION_BASED_STATISTICS =>'',
  :DATABASE =>'nr',
  #:DB_GENETIC_CODE =>'',
  :DESCRIPTIONS =>'1',
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
seq = ">gi|gatggghghg|embl|mnm,mnbjbkbkbjkbkkjb
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
  REGEX_NA=Regexp.new('[a|c|g|t|u]*')
  URL=URI.parse('http://www.ncbi.nlm.nih.gov/blast/Blast.cgi')
  #END----------CONSTANTS-------------------------------------------------
  
  #BEGIN--------ATTRIBUTES-------------------------------------------------
  
  # error_message - contains the information on errors which occurred 
  # during setting up and performing the Blast search
  attr_reader :error_message
  
  # seq_type - contains information of the type of the sequence provided
  # 'AA' - protein; 'NA' - nucleotide
  attr_reader :seq_type
  
  # header - contains information from a header of the sequence in FASTA format
  # could be changed if needed, after creation of the object 
  attr :header
  
  # rid - contains a unique RID code assigned during submission of the Blast search
  # used later on to retrive the data
  attr_reader :rid
  
  # done - contains information whether remote Blast search has been succesfully finished 
  # and data have be retrieved
  attr_reader :done
  
  # query - contains the String with an actual sequence, with removed not-allowed symbols,
  # wich has been submitted to Blast program
  attr_reader :query
  
  # @program - program used for Blast search
  @program
  
  # @error - contains information whether error has occured during setting up
  # and performing the Blast search
  @error=false
  
  @error_message=''
  @rid=''
   
  #BEGIN--------INITIALIZE-------------------------------------------------
  def initialize(sequence, params={}, seq_type_detect=false)
    
    # Read the parameters for the Blast program
    @query_params=
    {
      :PROGRAM => 'blastn'
    }.merge(params)
        
    @done=false
    
    # Reads QUERY sequence and checks its validity. If sequence (in FASTA or plain format) contains
    # not allowed symbols removes them and autodetects the sequence type.
    # The header of the sequence in FASTA format is stored in @header variable (otherwise @header is nil)
    seq_temp=sequence
    @header=seq_temp.slice!(/(^>.*\n)/).to_s
    @query=Bio::Sequence.auto(seq_temp)
    case @query.seq.class.to_s
      when "Bio::Sequence::AA" then
      @query=seq_temp.upcase.scan(REGEX_AA).to_s
      @seq_type='AA'
      when "Bio::Sequence::NA" then
       @query=seq_temp.downcase.scan(REGEX_NA).to_s
       @seq_type='NA'
    end
    if @query.to_s.length<1
      @error=true
      @error_message="NO VALID QUERY PROVIDED"
    end
    @query=@header + @query
    
    # Assigns BLAST PROGRAM - if sequence preprocessing succeeded
    # If no PROGRAM provided 'blastn' is being tried to be assigned
    # If PROGRAM does not apply to the QUERY type the error is returned
    if not @error then
    @program=@query_params[:PROGRAM]
    if not PROGRAMS.include?(@program)
      @error_message='UNKNOWN PROGRAM'
      @error=true
    elsif ((@seq_type=='NA' and not PROGRAMS[0..2].include?(@program)) or (@seq_type=='AA' and not PROGRAMS[3..4].include?(@program)))
      @error_message= "\'%s\' SEQUENCE PROVIDED FOR \'%s\' PROGRAM!" % [@seq_type, @program.upcase]
      @error=true
    end
  end
    # Reports error /for debugging only?/
    if @error
      puts "Problem with initialization: " + @error_message
    else do_blast
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
  private :do_blast
  #END----------DO_BLAST------------------------------------------------------

  #BEGIN--------FETCH_RESULTS-------------------------------------------------
  # Fetches the results of the commited blast search
  # Returns result of the search in XML format if successfully retrieved the data or false if error occurred
  # Accepts one optional argument:
  #                 true  - waits until BLAST search is finished and returns the the content of the output file
  #                         which can be used by Bio::Blast::Report.new() function to create a BLAST report
  #                 false - checks if BLAST search is finished if not reports that results are not yet available
  def fetch_results(wait=true)
    if @rid!=nil
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
puts h.rid
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
puts seq
puts h.query
#h.process_output