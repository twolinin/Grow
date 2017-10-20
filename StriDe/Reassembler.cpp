///-----------------------------------------------
// 
// Written by Jyun-Hong Lin
// 
//-----------------------------------------------
//
// PacBioScaffoldProcess.cpp - Scaffold of PacBio reads using FM-index
//

#include "Reassembler.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "ASQG.h"
#include "SequenceProcessFramework.h"
#include "ReassemblerProcess.h"
#include "BWTIntervalCache.h"
#include "ContigGraph.h"

#include "CGVisitors.h"

//
// Getopt
//
#define SUBPROGRAM "Reassembler"
static const char *CORRECT_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by JH.\n"
"\n"
"Copyright 2015 National Chung Cheng University\n";

static const char *CORRECT_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Correct PacBio reads via FM-index walk\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -p, --prefix=PREFIX              use PREFIX for the names of the index files (default: prefix of the input file)\n"
"      -o, --outfile=FILE               write the corrected reads to FILE (default: READSFILE.ec.fa)\n"
"      -t, --threads=NUM                use NUM threads for the computation (default: 1)\n"
"\nPacBio Scaffold parameters:\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 31, recommend: 31 (PacBioH), 17 (PacBioS).)\n"
"      -x, --kmer-threshold=N           Attempt to correct kmers that are seen less than N times. (default: 3)\n"
"      -y, --seed-kmer-threshold=N      Attempt to find kmers of seed that are seen large than N times. (default: 10)\n"
"      -d, --downward=N                 for each possible source, we consider N downward seeds as targets. (default: 1)\n"
"      -c, --collect=N                  for each possible source, we consider N downward seeds to collect reads. (default: 5)\n"
"      --first                          find overlap between two contigs\n"
"      --second                         find gap between two contigs\n"
"      --third                          use two TGS reads reassble contigs\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static std::string prefix;
    static std::string contigPrefix;
    static std::string readsFile;
    static std::string outFile;

    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;

    static int kmerLength = 15;
    static int kmerThreshold = 2;

    static int seedKmerThreshold = 2;
    static int collect = 5;
    static int tgsAvgLen = 3000;
	static int searchRange = 5000;
    
    static bool isFirst = false;
    static bool isSecond = false;
	static bool isThird = false;
}

static const char* shortopts = "p:t:o:a:k:x:y:d:c:v:P:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_FIRST, OPT_SECOND, OPT_THIRD};

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "threads",       required_argument, NULL, 't' },
    { "outfile",       required_argument, NULL, 'o' },
    { "prefix",        required_argument, NULL, 'p' },
    { "kmer-size",     required_argument, NULL, 'k' },
    { "kmer-threshold" ,required_argument, NULL, 'x' },
    { "seed-kmer-threshold"   ,required_argument, NULL, 'y' },
    { "downward"       ,required_argument, NULL, 'd' },
    { "collect"        ,required_argument, NULL, 'c' },
    { "contigPrefix"   ,required_argument, NULL, 'P' },
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { "first",         no_argument,       NULL, OPT_FIRST },
    { "second",        no_argument,       NULL, OPT_SECOND },
	{ "third",         no_argument,       NULL, OPT_THIRD },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int ReassemblerMain(int argc, char** argv)
{
    ReassemblerOptions(argc, argv);

    // Set the Reassembler parameters
    ReassemblerParameters ecParams;
    BWT *pBWT, *pRBWT;
    SampledSuffixArray* pSSA;

    BWT *cpBWT, *cpRBWT;
    SampledSuffixArray* cpSSA;
    
    ContigGraph* pGraph; 
    pGraph = CGUtil::loadFASTA(opt::readsFile);
    
    // Load indices
    #pragma omp parallel
    {
        #pragma omp single nowait
        {    //Initialization of large BWT takes some time, pass the disk to next job
            std::cout << std::endl << "Loading BWT:\t" << opt::prefix + BWT_EXT << "\n";
            pBWT = new BWT(opt::prefix + BWT_EXT, opt::sampleRate);
        }
        #pragma omp single nowait
        {
            std::cout << "Loading RBWT:\t" << opt::prefix + RBWT_EXT << "\n";
            pRBWT = new BWT(opt::prefix + RBWT_EXT, opt::sampleRate);
        }
        #pragma omp single nowait
        {
            std::cout << "Loading Sampled Suffix Array:\t" << opt::prefix + SAI_EXT << "\n";
            pSSA = new SampledSuffixArray(opt::prefix + SAI_EXT, SSA_FT_SAI);
        }
    }

    // Load indices
    #pragma omp parallel
    {
        #pragma omp single nowait
        {    //Initialization of large BWT takes some time, pass the disk to next job
            std::cout << std::endl << "Loading BWT:\t" << opt::contigPrefix + BWT_EXT << "\n";
            cpBWT = new BWT(opt::contigPrefix + BWT_EXT, opt::sampleRate);
        }
        #pragma omp single nowait
        {
            std::cout << "Loading RBWT:\t" << opt::contigPrefix + RBWT_EXT << "\n";
            cpRBWT = new BWT(opt::contigPrefix + RBWT_EXT, opt::sampleRate);
        }
        #pragma omp single nowait
        {
            std::cout << "Loading Sampled Suffix Array:\t" << opt::contigPrefix + SAI_EXT << "\n";
            cpSSA = new SampledSuffixArray(opt::contigPrefix + SAI_EXT, SSA_FT_SAI);
        }
    }
    
    std::string tgsFa    = opt::prefix + ".fa";
    std::string tgsFasta = opt::prefix + ".fasta";
	
    std::ifstream TGSin1( tgsFa.c_str() );
    std::ifstream TGSin2( tgsFasta.c_str() );	

	if(!TGSin1.is_open() && !TGSin2.is_open())
	{
		std::cout << "can't open TGS rawdata, please check prefix is fa or fasta.\n";
		std::cout << "...nothing to do\n";
		return 0;
	}
	
    std::string name;
    std::string data;
	
    size_t tgsReadNumber=0;
    size_t tgsTotalBase=0;
    
    while(TGSin1 && TGSin1.is_open())
    {
         TGSin1 >> name;
         TGSin1 >> data;
	 if( data.length() < 1000 ) continue;
         tgsTotalBase += data.length();
         tgsReadNumber++;
    }
	while(TGSin2 && TGSin2.is_open())
    {
         TGSin2 >> name;
         TGSin2 >> data;
	 if( data.length() < 1000 ) continue;
         tgsTotalBase += data.length();
         tgsReadNumber++;
    }
    opt::tgsAvgLen = tgsTotalBase/tgsReadNumber;
    
    BWTIndexSet indexSet;
    indexSet.pBWT = pBWT;
    indexSet.pRBWT = pRBWT;
    indexSet.pSSA = pSSA;
    ecParams.indices = indexSet;
    
    BWTIndexSet CindexSet;
    CindexSet.pBWT = cpBWT;
    CindexSet.pRBWT = cpRBWT;
    CindexSet.pSSA = cpSSA;
    ecParams.contigIndices = CindexSet;
    
    // Open outfiles and start a timer
    //std::ostream* pWriter = createWriter(opt::outFile);
    //std::ostream* pDiscardWriter = (!opt::discardFile.empty() ? createWriter(opt::discardFile) : NULL);
    Timer* pTimer = new Timer(PROGRAM_IDENT);

    ecParams.kmerLength = opt::kmerLength;
    ecParams.seedKmerThreshold = opt::seedKmerThreshold;
    ecParams.FMWKmerThreshold  = opt::kmerThreshold;
    ecParams.collectedSeeds    = opt::collect;
    ecParams.maxSeedInterval   = 500;
    ecParams.tgsAvgLen      = opt::tgsAvgLen*1.5;
	ecParams.searchRange    = opt::searchRange;
    ecParams.isFirst  = opt::isFirst;
    ecParams.isSecond = opt::isSecond;
    ecParams.isThird  = opt::isThird;
	
    std::cout << std::endl 
	      << "Reassemble reads for :    " << opt::readsFile              << std::endl
              << "number of threads    :    " << opt::numThreads             << std::endl
              << "large kmer size      :    " << ecParams.kmerLength         << std::endl 
              << "large kmer freq. cutoff : " << ecParams.seedKmerThreshold  << std::endl
              << "small kmer freq. cutoff : " << ecParams.FMWKmerThreshold   << std::endl
              << "TGS average length   :    " << opt::tgsAvgLen              << std::endl
              << "search range : "<< ecParams.searchRange << std::endl
	      << std::endl;
				
    // Setup post-processor
    ReassemblerPostProcess postProcessor(ecParams, pGraph);
    
    if(opt::numThreads <= 1)
    {
        // Serial mode
        ReassemblerProcess processor(ecParams);
        
        SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
                                                         ReassemblerResult,
                                                         ReassemblerProcess,
                                                         ReassemblerPostProcess>(opt::readsFile, &processor, &postProcessor);
    }
    else
    {
        // Parallel mode
        std::vector<ReassemblerProcess*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            ReassemblerProcess* pProcessor = new ReassemblerProcess(ecParams);
            processorVector.push_back(pProcessor);
        }

        SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
        ReassemblerResult,
        ReassemblerProcess,
        ReassemblerPostProcess>(opt::readsFile, processorVector, &postProcessor);

        for(int i = 0; i < opt::numThreads; ++i)
            delete processorVector[i];
    }

    postProcessor.filterErrorTGS();

    if(opt::isFirst)
    {
        postProcessor.buildGraphByOneTGS();
		
	    SGGraphStatsVisitor statsVisit;
        pGraph->visit(statsVisit);        
        pGraph->writeDot("1_Raw_Overlap_Graph.dot", 0);
        
        pGraph->Simplify();
        pGraph->writeDot("2_make_Overlap_Graph.dot", 0);
        
        SGFastaVisitor av("Reassembler_first_Overlap.fa");
        pGraph->visit(av);
    }
    else if(opt::isSecond)
    {
	    postProcessor.buildGraphByOneTGS();
		
	    SGGraphStatsVisitor statsVisit;
        pGraph->visit(statsVisit);        
        pGraph->writeDot("1_Raw_NonOverlap_Graph.dot", 0);
        
        pGraph->pacbioSimplify();
        pGraph->writeDot("2_make_NonOverlap_Graph.dot", 0);
        
        SGFastaVisitor av("Reassembler_second_NonOverlap.fa");
        pGraph->visit(av);
    }
    else if(opt::isThird)
    {
	    postProcessor.buildGraphByTwoTGS();
		
	    SGGraphStatsVisitor statsVisit;
        pGraph->visit(statsVisit);        
        pGraph->writeDot("1_Raw_NonOverlapByTwoTGS_Graph.dot", 0);
        
        pGraph->pacbioSimplify();
        pGraph->writeDot("2_make_NonOverlapByTwoTGS_Graph.dot", 0);
        
		//SGFastaVisitor av("Scaffold_third_NonOverlapByTwoPB.fa");
        SGFastaVisitor av("Reassembler_third_NonOverlapByTwoPB.fa");
        
		pGraph->visit(av);
    }
    
    delete pBWT;
    if(pRBWT != NULL)
    delete pRBWT;

    if(pSSA != NULL)
    delete pSSA;

    delete pTimer;

    return 0;
}


//
// Handle command line arguments
//
void ReassemblerOptions(int argc, char** argv)
{
    optind=1;    //reset getopt
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
        case 'p': arg >> opt::prefix; break;
        case 'o': arg >> opt::outFile; break;
        case 't': arg >> opt::numThreads; break;
        case 'k': arg >> opt::kmerLength; break;
        case 'x': arg >> opt::kmerThreshold; break;
        case '?': die = true; break;
        case 'v': opt::verbose++; break;
        case 'y': arg >> opt::seedKmerThreshold; break;
        case 'c': arg >> opt::collect; break;
        case 'P': arg >> opt::contigPrefix; break;
        case OPT_FIRST:  opt::isFirst  = true; break;
        case OPT_SECOND: opt::isSecond = true; break;
		case OPT_THIRD:  opt::isThird  = true; break;
        case OPT_HELP:
            std::cout << CORRECT_USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        case OPT_VERSION:
            std::cout << CORRECT_VERSION_MESSAGE;
            exit(EXIT_SUCCESS);
        
        }
    }

    if (argc - optind < 1)
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    }
    else if (argc - optind > 1)
    {
        
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if(opt::numThreads <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
        die = true;
    }

    if(opt::kmerLength <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid kmer length: " << opt::kmerLength << ", must be greater than zero\n";
        die = true;
    }

    if(opt::kmerThreshold <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid kmer threshold: " << opt::kmerThreshold << ", must be greater than zero\n";
        die = true;
    }
     
    if (die)
    {
        std::cout << "\n" << CORRECT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filenames
    opt::readsFile = argv[optind++];

    if(opt::prefix.empty())
    {
        opt::prefix = stripFilename(opt::readsFile);
    }

    // Set the Reassemble threshold
    if(opt::kmerThreshold <= 0)
    {
        std::cerr << "Invalid kmer support threshold: " << opt::kmerThreshold << "\n";
        exit(EXIT_FAILURE);
    }
    
    if( !opt::isFirst && !opt::isSecond && !opt::isThird)
    {
        std::cerr << "missing arguments --first or --second or --third"  << "\n";
        exit(EXIT_FAILURE);
    }
    
    if( (opt::isFirst && opt::isSecond) || (opt::isFirst && opt::isThird) || (opt::isSecond && opt::isThird) )
    {
        std::cerr << "Choice arguments --first or --second or --third"  << "\n";
        exit(EXIT_FAILURE);
    }

}
