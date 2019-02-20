#include "vcf_parser.h"
#include <Rcpp.h>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace {
    using namespace Rcpp;
    using namespace vcf;
    using namespace std;

    class Parser: public vcf::VCFParser {
        void handle_error(const vcf::ParserException& e) {
            Rf_warning(e.get_message().c_str());
        }
    public:
        using vcf::VCFParser::VCFParser;
    };
}

// [[Rcpp::export]]
List parse_vcf(const CharacterVector& filename, const CharacterVector& samples,
               const CharacterVector& bad_positions, const CharacterVector& available_positions,
               const IntegerVector& DP, const IntegerVector& GQ,
               const LogicalVector& ret_gmatrix, const CharacterVector& binary_prefix) {
    string name = (string)filename[0];
    bool gzipped = boost::algorithm::ends_with(name, ".vcf.gz");
    if(!gzipped && !boost::algorithm::ends_with(name, ".vcf")) {
        Rf_error("Filename must have extension .vcf or .vcf.gz");
    }
    VCFFilter filter(DP[0], GQ[0]);
}
