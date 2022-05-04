#pragma once

#include <functional>
#include "Types.h"

namespace maxwell {

class Source {
public:

	Source(FieldType&, std::function<double(const Position&)>);
	Source(FieldType&, std::function<double(const Position&)>, Direction&);

	FieldType& getFieldType() { return ft_; }
	std::function<double(const Position&)>& getFunction(){ return function_; }
	Direction& getDirection() { return d_; }

private:

	FieldType ft_;
	std::function<double(const Position&)> function_;
	Direction d_;

};
}