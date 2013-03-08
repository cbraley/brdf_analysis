#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <string>
#include <ostream>
#include <sstream>
#include <ios>

// Utility functions for string processing.
namespace StringUtils {

    // Check if a string ends with a particular suffix.
    bool endsWith(const std::string& to_check, const std::string& suffix);

    // Check if a string starts with a particular prefix.
    bool startsWith(const std::string& to_check, const std::string& prefix);

    // Parse an argument from arg_text, where the argument name is given by flag_text
    // The type T must be default constructible and have an overloaded stream insertion
    // operator.  Example usage:
    //  
    // bool parse_ok = true;
    // int num_days = StringUtils::parseArg(argv[i], "--num_days", &parse_ok);
    //
    // could be used to parse arguments where argv[i] is of the form
    // "--num_days=10".  In this case, parse_ok would be true.  However, for 
    // argv[i] of the form "--num_days10" or "--num_days 10" parse_ok
    // would be false and the return value num_days would be undefined.
    template<typename T>
    T parseArg(const std::string& arg_text, const std::string& flag_text,
               bool* parse_ok);

    // These are not public functions!  Unfortunately no way to hide their implementation 
    // since they are called by template functions...
    namespace Implementation {
        std::string extractAfterEquals(const std::string& text, const std::string& flag);
    }
}

// Implementation -------------------------------------------------------------
inline std::string StringUtils::Implementation::extractAfterEquals(
    const std::string& text,
    const std::string& flag)
{
    size_t loc = text.find(flag);
    if (loc == std::string::npos || loc + flag.size() > flag.size()) {
        return "";
    } else {
        return text.substr(loc + flag.size(), std::string::npos);
    }
}

template<typename T>
inline T StringUtils::parseArg(
    const std::string& arg_text,
    const std::string& flag_text,
    bool* parse_ok) 
{
    std::string arg_val_txt;
    arg_val_txt = StringUtils::Implementation::extractAfterEquals(
        arg_text, flag_text + "=");
    T ret;
    if (arg_val_txt.empty()) {
        *parse_ok = false;
        return ret;
    } else {
        std::stringstream ss(arg_val_txt);
        ss << std::boolalpha;
        ss >> ret;
        *parse_ok = !ss.fail();
        return ret;
    }
}



#endif //STRING_UTILS_H
