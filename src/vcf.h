#ifndef SRC_VCF_H
#define SRC_VCF_H

#include <string>
#include <vector>
#include <istream>
#include <unordered_map>
#include <unordered_set>
#include <functional>

#include <boost/functional/hash.hpp>

namespace vcf {
    class Variant;
}

namespace vcf {

    class ParserException {
        std::string msg;
    public:
        explicit ParserException(const std::string& message);

        std::string get_message();
    };

    class Chromosome {
        const int chrX = 23;
        const int chrY = 24;
        int chr;

        bool parse(std::string);

    public:
        Chromosome(const std::string&);
        explicit operator std::string() const;
        int num() const;

        friend bool operator==(const Chromosome& chr, const Chromosome& other);
    };

    class Position {
        Chromosome chr;
        int pos;
    public:
        Position(Chromosome chr, int pos);
        Position(const Position& other);
        explicit operator std::string() const;
        Chromosome chromosome() const;
        int position() const;

        friend bool operator==(const Position& pos, const Position& other);
        friend size_t hash_value(const Position& pos);
    };


    class Variant {
        Position pos;
        std::string ref;
        std::string alt;
    public:
        Variant(Position pos, std::string& ref, std::string& alt);
        explicit operator std::string() const;
        Position position() const;
        std::string reference() const;
        std::string alternative() const;

        friend bool operator==(const Variant& variant, const Variant& other);
    };
}

namespace std {
    template <>
    struct hash<vcf::Position> {
        size_t operator()(const vcf::Position& pos) const {
            return hash_value(pos);
        }
    };

    template <>
    struct hash<vcf::Variant> {
        size_t operator()(const vcf::Variant& var) const {
            size_t seed = 0;
            boost::hash_combine(seed, var.position());
            boost::hash_combine(seed, var.reference());
            boost::hash_combine(seed, var.alternative());
            return seed;
        }
    };
}

namespace vcf {
    class Range {
        Chromosome chr;
        int from;
        int to;
    public:
        Range(Chromosome chr, int from, int to);
        bool includes(const Position& p) const;
    };

    class VCFFilter {
        std::unordered_set<Variant> available_variants;
        std::unordered_set<Variant> bad_variants;
        std::unordered_set<std::string> available_samples;
        std::vector<Range> ranges;
    public:
        VCFFilter();
        bool apply(const Variant& v) const;
    };

    class VCFParser {
        std::vector<std::string> fields;
        std::vector<Variant> variants_found;
        std::vector<std::string> samples_found;
        std::vector<std::vector<int>> gmatrix;
        std::unordered_map<std::string, int> fields_positions;

    public:
        explicit VCFParser(std::istream stream, VCFFilter& filter);
};

}

#endif //SRC_VCF_H
