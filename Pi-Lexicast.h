#ifndef PI_LEXICAST_H_INCLUDED
#define PI_LEXICAST_H_INCLUDED
#include<string>
#include<sstream>
#include<cstdlib>
#include<stdexcept>
#include<limits>
#include<mutex>

namespace Pitty
{

struct lexicast_static_istringstream_
{
    static std::istringstream& get()
    {
        static std::istringstream iss;
        return iss;
    }

    //static std::mutex my_mutex;
};
// generic string -> T converter
template <typename T>
T lexicast(const std::string& str)
{
    T var;
    // there is a unique static istringstream object : it is not thread safe
    static std::istringstream& iss=lexicast_static_istringstream_::get();

    iss.str(str);
    iss >> var;
    // deal with any error bits that may have been set on the stream
    return var;
}



template <>
double lexicast<double>(const std::string& str)
{
    return atof(str.c_str());
}


template <>
int lexicast<int>(const std::string& str)
{
    return atoi(str.c_str());
}
template <>
unsigned lexicast<unsigned>(std::string const & str)
 {
    long result = atol(str.c_str());
    if (result > std::numeric_limits<unsigned>::max() || result <0)
    {
        throw std::out_of_range("lexicast<unsigned>");
    }
    return unsigned(result);
}








}// end namespace


#endif // PI_LEXICAST_H_INCLUDED
