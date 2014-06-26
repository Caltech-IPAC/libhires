#ifndef HIRES_EXCEPTION_HXX
#define HIRES_EXCEPTION_HXX

namespace hires
{
  class Exception: public std::runtime_error
  {
  public:
    Exception(const std::string &s): std::runtime_error(s) {}
  };
}


#endif
