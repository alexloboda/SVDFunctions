#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <stdexcept>
#include <string>

namespace mvn {

class MatrixException : public std::logic_error {
public:
    explicit MatrixException(const std::string& message) : std::logic_error(message) {}
};

class ClusteringException : public std::invalid_argument {
public:
    explicit ClusteringException(const std::string& message) : std::invalid_argument(message) {}
};

class SamplingException : public std::logic_error {
public:
    explicit SamplingException(const std::string& message) : std::logic_error(message) {}
};

} // namespace mvn

#endif // EXCEPTIONS_H
