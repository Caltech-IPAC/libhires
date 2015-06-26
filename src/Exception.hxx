#pragma once

#include <stdexcept>

namespace hires
{
class Exception : public std::runtime_error
{
public:
  Exception (const std::string &s) : std::runtime_error (s) {}
};
}
