require 'ncbi_blast.class.rb'

#Initialization hash - comment or remove parameters which are not going to be modified
params=
{
  :ALIGNMENTS =>'10',
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
  #:FILTER =>'L',
  #:FORMAT_ENTREZ_QUERY =>'',
  #:FORMAT_OBJECT =>'',
  #:FORMAT_TYPE =>'',
  #:GAPCOSTS =>'',
  #:GENETIC_CODE =>'',
  #:HITLIST_SIZE =>'',
  #:I_THRESH =>'',
  #:LAYOUT =>'',
  #:LCASE_MASK =>'',
  # :MEGABLAST =>'yes',
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
seq = ">gi|gatggghghg|embl|mnm,mnbjbkbkbjkbkkjb
ATGcagctctttgtccgcgcccaggagctacacaccttcgaggtgaccggccaggaaacg
gtcgcccagatcaaggctcatgtagcctcactggagggcattgccccggaagatcaagtc
gtgctcctggcaggcgcgcccctggaggatgaggccactctgggccagtgcggggtggag
gccctgactaccctggaagtagcaggccgcatgcttggaggtaaagtccatggttccctg
gcccgtgctggaaaagtgagaggtcagactcctaaggtggccaaacaggagaagaagaag
aagaagacaggtcgggctaagcggcggatgcagtacaaccggcgctttgtcaacgttgtg
cccacctttggcaagaagaagggccccaatgccaactcttaa"



# Putting BLAST query onto the host. Making the PUT BLAST Command.
# Constructor also waits till the query is completed on the host. (see definition of constructor)
tries = 3
begin
  h=NCBI_BLAST.new(seq, params)
rescue StandardError,Timeout::Error,OpenURI::HTTPError => e
  tries -= 1
  puts "Error message: " + e.message + " . " + tries.to_s + " tries left."
  if tries >= 0
    sleep 60
    retry
  end
end

# If tries >= 0, then the constructor method succeeded eventually.
if tries >= 0
  puts "Request ID: " + h.rid

  # Fetching results from the host. Build Report object from retrieved data.
  tries = 6
  begin
    rep=Bio::Blast::Report.new(h.fetch_results)
  rescue StandardError,Timeout::Error,OpenURI::HTTPError => e
    tries -= 1
    puts "Error message: " + e.message + " . " + tries.to_s + " tries left."
    if tries >= 0 && e.message != "INVALID RID"
      sleep 60
      retry
    end
  end
end

# Process results
# See BioRuby documentation for details
if not rep.nil?
  rep.each do |h|
    puts h.definition + "\n\r"
  end
end
