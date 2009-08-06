#---------------------------------------------------
require 'rubygems'
require 'uri'
require 'logger'
require 'bio'
require 'open-uri'
require 'net/http'

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

  # answer - contains the result of the query
  attr_reader :answer

  # query - contains the String with an actual sequence, with removed not-allowed symbols,
  # wich has been submitted to Blast program
  attr_reader :query

  # @program - program used for Blast search
  @program

  # @logger - an instance of the ruby standart library logger class used for logging.
  @logger

  # @error - contains information whether error has occured during setting up
  # and performing the Blast search
  @error=false

  @error_message=''
  @rid=''

  #BEGIN--------INITIALIZE-------------------------------------------------
  def initialize(sequence, params={},logger=Logger.new(STDOUT), seq_type_detect=false)

    @logger = logger

    # Read the parameters for the Blast program
    @query_params=
    {
      :PROGRAM => 'blastn'
      }.merge(params)


      # Reads QUERY sequence and checks its validity. If sequence (in FASTA or plain format) contains
      # not allowed symbols removes them and autodetects the sequence type.
      # The header of the sequence in FASTA format is stored in @header variable (otherwise @header is nil)
      seq_temp=sequence
      @header=seq_temp.slice!(/(^>.*\n)/).to_s
      @query=Bio::Sequence.new(seq_temp)
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
        @logger.error "Problem with initialization: " + @error_message
        raise StandardError.new("Problem with initialization: " + @error_message)
      else
        # Do BLAST search and wait for results.
        # Maybe we should move this code to a new function.
        do_blast
        if @rid.length > 0
          res='http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=' + @rid
          while status = open(res).read.scan(/Status=(.*?)$/).to_s=='WAITING'
            @logger.debug("Status=WAITING")
            sleep(3)
          end  
          case status
          when 'FAILED' then 
            @logger.error "Fetching of data failed. RID: " + @rid.to_s
            raise StandardError.new("Fetching of the data failed")
          when 'UNKNOWN' then
            @logger.error "RID UNKNOWN OR EXPIRED. RID: " + @rid.to_s
            raise StandardError.new("RID UNKNOWN OR EXPIRED")
          when 'READY' then
            @logger.debug("Status=READY." + " RID: " + @rid.to_s )
          end
        else
          @logger.error "BLAST SEARCH DID NOT INITIATED. RID == nil."
          raise StandardError.new("BLAST SEARCH DID NOT INITIATED. RID == nil.") 
        end
      end
    end
    #END----------INITIALIZE-------------------------------------------------

    #BEGIN--------DO_BLAST-------------------------------------------------
    # Returns String value representing assigned RID
    # In case of the error returns nil

    def do_blast
      sub=@query_params.merge({'CMD'=>'Put', 'QUERY'=>@query})
      begin
        http = Net::HTTP::post_form(URL, sub)
      
      # Catch a whole bunch of exceptions, that might raised by the post_form method. 
      # Simplifies caller's error handling, we don't need detailed information about the HTTP error.
      rescue Timeout::Error, Errno::EINVAL, Errno::ECONNRESET, EOFError, Net::HTTPBadResponse, Net::HTTPHeaderSyntaxError, Net::ProtocolError => e
        @logger.error "Putting BLAST search query onto host failed. Put Command failed."
        raise StandardError.new("Putting BLAST search query onto host failed. Description: #{e.message}")
      end
      @rid=http.body.scan(/RID = (.*?)$/).to_s
    end
    private :do_blast
    #END----------DO_BLAST------------------------------------------------------

    #BEGIN--------FETCH_RESULTS-------------------------------------------------
    # Fetches the results of the commited blast search
    # Returns result of the search in XML format if successfully retrieved the data.
    # Exception handling must done by the caller.
    def fetch_results
      if @rid.length > 0
        res='http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID=' + rid.to_s
        @answer=open(res).read
        return @answer
      else 
        @logger.error "BLAST SEARCH DID NOT INITIATED. RID == nil."
        raise StandardError.new("INVALID RID")
      end
    end
    #END----------FETCH_RESULTS------------------------------------------------
  end
  #------END of NCBI_BLAST CLASS DEFINITION------------------------------------