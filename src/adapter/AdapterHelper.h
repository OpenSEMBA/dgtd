#pragma once
#include <string>

inline void checkIfThrows(bool condition, const std::string& msg)
{
	if (!condition) {
		throw std::runtime_error(msg.c_str());
	}
}
