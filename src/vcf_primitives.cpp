#include "vcf_primitives.h"

#include <string>
#include <algorithm>

namespace {
    using std::string;
    using std::to_string;
    using std::transform;
    using std::invalid_argument;
}

namespace vcf {
    Position::Position(Chromosome chr, int pos) : pos(pos), chr(chr) {}

    Position::Position(const vcf::Position& other)
            : chr(other.chr), pos(other.pos) {}

    bool operator==(const Position& pos, const Position& other) {
        return pos.chr == other.chr && pos.pos == other.pos;
    }

    Position::operator std::string() const {
        return (string)chr + ":" + to_string(pos);
    }

    Chromosome Position::chromosome() const {
        return chr;
    }

    int Position::position() const {
        return pos;
    }

    size_t hash_value(const Position& pos) {
        size_t seed = 0;
        boost::hash_combine(seed, pos.chromosome().num());
        boost::hash_combine(seed, pos.position());
        return seed;
    }

    Variant::Variant(Position pos, const string& ref, const string& alt)
            : pos(pos), ref(ref), alt(alt) {}

    bool operator==(const Variant& var, const Variant& other) {
        return var.pos == other.pos && var.ref == other.ref && var.alt == other.alt;
    }

    Variant::operator std::string() const {
        return (string)pos + "\t" + ref + "\t" + alt;
    }

    Position Variant::position() const {
        return pos;
    }

    std::string Variant::reference() const {
        return ref;
    }

    std::string Variant::alternative() const {
        return alt;
    }

    bool Chromosome::parse(std::string str) {
        if (!std::all_of(str.begin(), str.end(), ::isdigit)) {
            transform(str.begin(), str.end(), str.begin(), ::tolower);
            if (str.substr(0, 3) != "chr") {
                return false;
            }
            str = str.substr(3, str.length());
        }

        if (str == "X") {
            chr = chrX;
        }
        if (str == "Y") {
            chr = chrY;
        }

        try {
            chr = stoi(str);
            if (chr < 1 || chr > 22) {
                return false;
            }
        } catch (invalid_argument& e) {
            return false;
        }
        return true;
    }

    Chromosome::Chromosome(const std::string& str) {
        if (!parse(str)) {
            throw ParserException(R"(Parser error: expected "chr##", found )" + str);
        }
    }

    bool operator==(const Chromosome& chr, const Chromosome& other) {
        return chr.chr == other.chr;
    }

    Chromosome::operator std::string() const {
        string str_rep = "chr";
        if (chr == chrX) {
            str_rep += "X";
        } else if ( chr == chrY) {
            str_rep += "Y";
        } else {
            str_rep += to_string(chr);
        }
        return str_rep;
    }

    int Chromosome::num() const {
        return chr;
    }

    ParserException::ParserException(std::string message) :msg(message) {}

    std::string ParserException::get_message() const {
        return msg;
    }

    ParserException::ParserException(std::string message, int line) {
        msg = "Error in line " + to_string(line) + ": " + message;
    }

    Range::Range(Chromosome chr, int from, int to)
        :chr(chr), from(from), to(to) {}

    bool Range::includes(const Position& p) const {
        if (!(chr == p.chromosome())) {
            return false;
        }
        return p.position() >= from && p.position() < to;
    }

    Allele::Allele(AlleleType type, int DP, int GQ) :type(type), depth(DP), quality(GQ){}

    int Allele::DP() const {
        return depth;
    }

    int Allele::GQ() const {
        return quality;
    }

    AlleleType Allele::alleleType() const {
        return type;
    }
}
