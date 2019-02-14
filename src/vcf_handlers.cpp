#include "vcf_handlers.h"

namespace {
    using std::vector;
}

namespace vcf {
    VariantsHandler::VariantsHandler(const std::vector<std::string>& samples) :samples(samples){}
}
