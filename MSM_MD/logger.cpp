#include "logger.h"

using namespace MSM_MD_NS;

FILE *Logger::screen = stdout;
FILE *Logger::logfile = nullptr;


void Logger::log(std::string str, bool to_screen, bool to_log, const char* newline) {
	str += newline;

	if (to_screen) {
		fputs(str.c_str(), screen);
	}
	if (to_log && logfile != nullptr)
	{
		fputs(str.c_str(), logfile);
	}
}

void Logger::error(std::string str, bool to_screen, bool to_log, const char* newline)
{
	std::string err_msg = "[ERROR] ";
	err_msg += str;

	Logger::log(err_msg, to_screen, to_log, newline);

	fclose(Logger::logfile);

	throw std::runtime_error(str);
}

void Logger::warning(std::string str, bool to_screen, bool to_log, const char* newline)
{
	std::string warn_msg = "[WARNING] ";
	warn_msg += str;

	Logger::log(warn_msg, to_screen, to_log, newline);
}

void Logger::log(char* str, bool to_screen, bool to_log, const char* newline)
{
	std::string chr_str = std::string(str);

	Logger::log(chr_str, to_screen, to_log, newline);
}

void Logger::error(char* str, bool to_screen, bool to_log, const char* newline)
{
	std::string chr_str = std::string(str);

	Logger::error(chr_str, to_screen, to_log, newline);
}

void Logger::warning(char* str, bool to_screen, bool to_log, const char* newline)
{
	std::string chr_str = std::string(str);

	Logger::warning(chr_str, to_screen, to_log, newline);
}

void Logger::log(char c, bool to_screen, bool to_log, const char* newline)
{
	std::string chr_str = "";

	Logger::log(chr_str + c, to_screen, to_log, newline);
}

void Logger::error(char c, bool to_screen, bool to_log, const char* newline)
{
	std::string chr_str = "";

	Logger::error(chr_str + c, to_screen, to_log, newline);
}

void Logger::warning(char c, bool to_screen, bool to_log, const char* newline)
{
	std::string chr_str = "";

	Logger::warning(chr_str + c, to_screen, to_log, newline);
}

void Logger::change_file(std::string filepath)
{

	logfile = fopen(filepath.c_str(), "w");
}