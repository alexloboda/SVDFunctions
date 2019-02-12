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
        static const int chrX = 23;
        static const int chrY = 24;

        int chr;

        bool parse(std::string);

    public:
        explicit Chromosome(const std::string&);
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

    class Range {
        Chromosome chr;
        int from;
        int to;
    public:
        Range(Chromosome chr, int from, int to);
        bool includes(const Position& p) const;
    };

    enum AlleleType {HOM, HET, ALT, MISSING};

    class Allele {
        int depth;
        int quality;
        AlleleType type;
    public:
        Allele(AlleleType type, int DP, int GQ);
        int DP() const;
        int GQ() const;
        AlleleType alleleType() const;
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

#endif //SRC_VCF_H
