#include "utils.h"

using namespace MSM_MD_NS;

bool utils::isInteger(std::string str) {
	if (str.empty())
	{
		return false;
	}

	try {
		long str_int = std::stol(str);
		return true;
	}
	catch (std::exception &e){
		return false;
	}
}

bool utils::isInteger(char* str)
{
	return utils::isInteger(std::string(str));
}

bool utils::isNumeric(std::string str) {
	if (str.empty())
	{
		return false;
	}

	try {
		double str_double = std::stod(str);
		return true;
	}
	catch (std::exception& e) {
		return false;
	}
}

bool utils::isNumeric(char* str)
{
	return utils::isNumeric(std::string(str));
}

long utils::toInteger(std::string str) {
	bool is_int = utils::isInteger(str);
	
	if (!is_int) {
		Logger::error("string \"" + str + "\" cannot be converted to integer");
	}
	return std::stol(str);
}

long utils::toInteger(char* str)
{
	return utils::toInteger(std::string(str));
}

double utils::toFloat(std::string str) {
	bool is_float = utils::isNumeric(str);

	if (!is_float) {
		Logger::error("string \"" + str + "\" cannot be converted to float");
	}
	return std::stod(str);
}

double utils::toFloat(char* str)
{
	return utils::toFloat(std::string(str));
}

std::vector<std::string> utils::split(std::string str) {
	std::vector<std::string> split_str = {};

	std::string curr_str;

	for (int i = 0; i < str.size(); i++) {
		if (!isspace(str[i]))
		{
			curr_str.push_back(str[i]);
		}
		else {
			if (!curr_str.empty())
			{
				split_str.push_back(curr_str);
			}
			curr_str = "";
		}
	}

	split_str.push_back(curr_str);

	return split_str;
}

std::vector<std::string> utils::split(std::string str, std::string sep) {
	std::vector<std::string> split_str = {};

	int f_ind = str.find(sep);
	std::string sub = str.substr(0, -1);

	int count = 0;

	while (f_ind >= 0)
	{
		split_str.push_back(sub.substr(0, f_ind));

		sub = sub.substr(f_ind + sep.size(), -1);

		f_ind = sub.find(sep);
	}

	split_str.push_back(sub.substr(0, f_ind));

	return split_str;
}

std::string utils::strip(std::string str) {
	int low = 0;

	for (int i = 0; i < str.size(); i++)
	{
		if (!isspace(str[i])) {
			low = i;
			break;
		}
	}

	int high = str.size() - 1;

	for (int i = str.size() - 1; i >= low; i--)
	{
		if (!isspace(str[i])) {
			high = i + 1;
			break;
		}
	}
	
	return str.substr(low, high - low);
}