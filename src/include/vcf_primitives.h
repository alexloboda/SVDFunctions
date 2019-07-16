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
        explicit ParserException(std::string message);
        ParserException(std::string message, int line);

        std::string get_message() const;
    };

    class Chromosome {
        static const int chrX = 23;
        static const int chrY = 24;

        int chr;

        bool parse(std::string);

    public:
        explicit Chromosome(int num);
        explicit Chromosome(const std::string&);
        Chromosome(const Chromosome& other) = default;
        Chromosome(Chromosome&& other) = default;
        Chromosome& operator=(const Chromosome& other) = default;
        Chromosome& operator=(Chromosome&& other) noexcept = default;
        explicit operator std::string() const;
        int num() const;

        friend bool operator==(const Chromosome& chr, const Chromosome& other);
        friend bool operator!=(const Chromosome& chr, const Chromosome& other);
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
        Position operator+(int shift);

        static Position parse_position(const std::string& str);

        friend bool operator==(const Position& pos, const Position& other);
        friend size_t hash_value(const Position& pos);
    };


    class Variant {
        Position pos;
        std::string ref;
        std::string alt;
    public:
        Variant(Position pos, const std::string& ref, const std::string& alt);
        explicit operator std::string() const;
        Position position() const;
        std::string reference() const;
        std::string alternative() const;

        static std::vector<Variant> parseVariants(const std::string& s);

        friend bool operator==(const Variant& variant, const Variant& other);
    };

    class Range {
        Chromosome chr;
        int from;
        int to;
    public:
        Range(Chromosome chr, int from, int to);
        Range(const Range& other) = default;
        Range(Range&& other) = default;
        Range& operator=(const Range& other) noexcept = default;
        Range& operator=(Range&& other) noexcept = default;
        bool includes(const Position& p) const;
        bool operator<(const Range& other) const;
        bool operator<(const Position& pos) const;
        Position begin() const;
        Position end() const;
        explicit operator std::string() const;

        static Range parseRange(const std::string& s);
    };

    class RangeSet {
        std::unordered_map<int, std::set<Range>> rs;
    public:
        RangeSet() = default;
        bool includes(const Position& position);
        void insert(Range r);
        bool empty() const;
    };

    enum AlleleType {HOMREF = 0, HET = 1, HOM = 2, MISSING = 3};
    int to_int(AlleleType allele);

    class Allele {
        unsigned depth;
        unsigned quality;
        AlleleType type;
    public:
        Allele(AlleleType type, unsigned DP, unsigned GQ);

        unsigned DP() const;
        unsigned GQ() const;
        AlleleType alleleType() const;
    };

    struct BinaryAllele {
        uint16_t DP;
        uint8_t GQ;
        uint8_t allele;
        friend std::ostream& operator<<(std::ostream& os, const BinaryAllele& dt);
        friend std::istream& operator>>(std::istream& in, BinaryAllele& obj);
        static BinaryAllele fromAllele(const Allele& allele);
        static Allele toAllele(const BinaryAllele& allele);
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
