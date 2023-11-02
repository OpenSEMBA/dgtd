#pragma once
#include <string>

inline void checkIfThrows(bool condition, const std::string& msg)
{
	if (!condition) {
		throw std::exception(msg.c_str());
	}
}