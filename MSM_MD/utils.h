#ifndef UTILS
#define UTILS

#include <string>
#include <vector>
#include "logger.h"


namespace MSM_MD_NS {
	namespace utils {
		bool isNumeric(std::string str); // Checks if a string can be converted to a float
		bool isInteger(std::string str); // Checks if a string can be converted to an integer
		double toFloat(std::string str); // Convert a string to a float
		long toInteger(std::string str); // Convert a string to an int

		bool isNumeric(char* str); // Checks if a string can be converted to a float
		bool isInteger(char* str); // Checks if a string can be converted to an integer
		double toFloat(char* str); // Convert a string to a float
		long toInteger(char* str); // Convert a string to an int

		std::vector<std::string> split(std::string str);
		std::vector<std::string> split(std::string str, std::string sep);

		std::string strip(std::string str);
	}
}
#endif