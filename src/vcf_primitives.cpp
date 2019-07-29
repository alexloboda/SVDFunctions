#include "include/vcf_primitives.h"

#include <string>
#include <algorithm>

namespace {
    using std::string;
    using std::to_string;
    using std::transform;
    using std::invalid_argument;
    using std::vector;

    using std::size_t;
}

namespace vcf {
    Position::Position(Chromosome chr, int pos) : chr(chr), pos(pos) {}

    Position::Position(const vcf::Position& other)
            : chr(other.chr), pos(other.pos) {}

    bool operator==(const Position& pos, const Position& other) {
        return pos.chr == other.chr && pos.pos == other.pos;
    }

    int to_int(vcf::AlleleType type) {
        switch(type) {
            case vcf::HOMREF:
                return 0;
            case vcf::HET:
                return 1;
            case vcf::HOM:
                return 2;
            case vcf::MISSING:
                return 3;
            default:
                throw std::logic_error("[internal] unknown allele type");
        }
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

    Position Position::parse_position(const std::string& str) {
        size_t pos = std::find(str.begin(), str.end(), ':') - str.begin();
        if (pos >= str.length() - 1 || pos == 0) {
            throw ParserException("Position must be in format chr#:# but " + str + " given");
        }
        Chromosome chr(str.substr(0, pos));
        try {
            size_t final;
            int position = std::stoi(str.substr(pos + 1, str.length() - pos), &final);
            if (final != str.length() - pos - 1) {
                throw std::invalid_argument("");
            }
            return {chr, position};
        } catch (...) {
            throw ParserException("Position must be in format chr#:# but \"" + str + "\" given");
        }
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

    vector<Variant> Variant::parseVariants(const std::string& s) {
        std::istringstream iss(s);
        vector<Variant> ret;
        if (s.empty()) {
            return ret;
        }
        string token;
        iss >> token;
        Position pos = Position::parse_position(token);
        string ref;
        iss >> ref;
        string alt;
        while(std::getline(iss >> std::ws, alt, ',')) {
            ret.emplace_back(pos, ref, alt);
        }
        return ret;
    }

    bool Chromosome::parse(std::string str) {
        transform(str.begin(), str.end(), str.begin(), ::tolower);
        if (str.substr(0, 3) == "chr") {
            str = str.substr(3, str.length());
        }

        if (str == "x") {
            chr = chrX;
            return true;
        }
        if (str == "y") {
            chr = chrY;
            return true;
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

    bool operator!=(const Chromosome& chr, const Chromosome& other) {
        return chr.chr != other.chr;
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

    Chromosome::Chromosome(int num) :chr(num){}

    ParserException::ParserException(std::string message) :msg(std::move(message)) {}

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

    Position Range::begin() const {
        return {chr, from};
    }

    Position Range::end() const {
        return {chr, to};
    }

    Position Position::operator+(int shift) {
        return {chr, pos + shift};
    }

    Range Range::parseRange(const std::string& s) {
        std::istringstream iss(s);
        string chr, start, end;
        int startpos = 0;
        int endpos = 0;
        try {
            iss >> chr >> start >> end;
            if (iss.fail()) {
                throw ParserException("Malformed range: " + s);
            }
            startpos = stoi(start);
            endpos = stoi(end) + 1;
        } catch(...) {
            throw ParserException("Malformed range: " + s);
        }
        Chromosome chromosome = Chromosome(chr);
        return {chromosome, startpos, endpos};
    }

    bool Range::operator<(const Range& other) const {
        return chr.num() < other.chr.num() || (chr == other.chr && to < other.to);
    }

    bool Range::operator<(const Position& pos) const {
        return chr.num() < pos.chromosome().num() || (chr == pos.chromosome() && to < pos.position());
    }

    Range::operator std::string() const {
        return "chr" + std::to_string(chr.num()) + " " + std::to_string(from) + " " + std::to_string(to - 1);
    }

    Allele::Allele(AlleleType type, unsigned DP, unsigned GQ) :depth(DP), quality(GQ), type(type) {}

    unsigned Allele::DP() const {
        return depth;
    }

    unsigned Allele::GQ() const {
        return quality;
    }

    AlleleType Allele::alleleType() const {
        return type;
    }

    std::ostream& operator<<(std::ostream& out, const BinaryAllele& allele) {
        out.write(reinterpret_cast<const char*>(&allele), sizeof(BinaryAllele));
        return out;
    }

    BinaryAllele BinaryAllele::fromAllele(const Allele& allele) {
        BinaryAllele binary{};
        unsigned DP = std::min(allele.DP(), 65535U);
        unsigned GQ = std::min(allele.GQ(), 255U);
        binary.DP = (uint16_t)DP;
        binary.GQ = (uint8_t)GQ;
        binary.allele = allele.alleleType();
        return binary;
    }

    std::istream& operator>>(std::istream& in, BinaryAllele& obj) {
        in.read((char*)&obj, sizeof(BinaryAllele));
        return in;
    }

    Allele BinaryAllele::toAllele(const BinaryAllele& allele) {
        auto type = (AlleleType)allele.allele;
        return {type, allele.DP, allele.GQ};
    }

    bool RangeSet::includes(const Position& pos) {
        auto& set = rs[pos.chromosome().num()];
        // Workaround, c++11
        auto it = set.lower_bound(Range(pos.chromosome(), -1, pos.position()));
        return it != set.end() && it->includes(pos);
    }

    void RangeSet::insert(Range r) {
        auto& ranges = rs[r.begin().chromosome().num()];
        auto pos = ranges.lower_bound(r);
        if (pos != ranges.end() && pos->includes(r.end() + (-1))) {
            int from = std::min(r.begin().position(), pos->begin().position());
            int to = std::max(r.end().position(), pos->end().position());
            r = Range(r.begin().chromosome(), from, to);
            ranges.erase(*pos);
        }
        ranges.insert(r);
    }

    bool RangeSet::empty() const {
        return rs.empty();
    }
}
