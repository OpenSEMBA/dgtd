#include "Sources.h"

namespace maxwell {

Source::Source(FieldType& ft, std::function<double(const Position&)> f) :
ft_(ft),
function_(f)	
{
}

Source::Source(FieldType& ft, std::function<double(const Position&)> f, Direction& d) :
ft_(ft),
function_(f),
d_(d)
{
}

}