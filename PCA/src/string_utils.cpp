#include "string_utils.h"

bool StringUtils::endsWith(const std::string& to_check, const std::string& suffix) {
    if (to_check.size() < suffix.size()) {
        return false;
    } else {
        return to_check.substr(to_check.size() - suffix.size(), suffix.size()) == suffix;
    }
}

bool StringUtils::startsWith(const std::string& to_check, const std::string& prefix) {
    return to_check.find(prefix) == 0;
}

