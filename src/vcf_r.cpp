#include "vcf_parser.h"
#include <Rcpp.h>
#include <boost/algorithm/string/predicate.hpp>
#include <iostream>
#include <fstream>
#include "zstr/zstr.hpp"
#include "zstr/strict_fstream.hpp"

namespace {
    using namespace Rcpp;
    using namespace vcf;
    using namespace std;
    using boost::algorithm::ends_with;

    class Parser: public VCFParser {
        void handle_error(const vcf::ParserException& e) {
            Rf_warning(e.get_message().c_str());
        }
    public:
        using VCFParser::VCFParser;
    };

    class RGenotypeMatrixHandler: public GenotypeMatrixHandler {
    public:
        using GenotypeMatrixHandler::GenotypeMatrixHandler;
        List result() {
            IntegerMatrix res(gmatrix.size(), samples.size());
            for (int i = 0; i < gmatrix.size(); i++) {
                for (int j = 0; j < samples.size(); j++) {
                    int val = gmatrix[i][j];
                    if (val == vcf::MISSING) {
                        val = NA_INTEGER;
                    }
                    res[i * samples.size() + j] = val;
                }
            }
            List ret;
            ret["gmatrix"] = res;
            return ret;
        }
    };
}

//' @export
// [[Rcpp::export]]
List parse_vcf(const CharacterVector& filename, const CharacterVector& samples,
               const CharacterVector& bad_positions, const CharacterVector& available_positions,
               const IntegerVector& DP, const IntegerVector& GQ,
               const LogicalVector& ret_gmatrix, const CharacterVector& binary_prefix) {
    const char* name = filename[0];
    unique_ptr<std::istream> in(new zstr::ifstream(name));
    VCFFilter filter(DP[0], GQ[0]);
    Parser parser(*in, filter);
    parser.parse_header();
    auto ss = parser.sample_names();
    auto hdlr = make_shared<RGenotypeMatrixHandler>(ss);
    parser.register_handler(hdlr);
    parser.parse_genotypes();
    return hdlr->result();
}
