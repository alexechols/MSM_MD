#ifndef LOGGER
#define LOGGER

#include <string>
#include <iostream>
#include <stdexcept>

namespace MSM_MD_NS {
	static class Logger {
	public: 
		static FILE *logfile;
		static FILE *screen;

		static void error(std::string str, bool to_screen = true, bool to_log = true, const char *newline="\n");
		static void warning(std::string str, bool to_screen = true, bool to_log = true, const char* newline = "\n");
		static void log(std::string str, bool to_screen = true, bool to_log = true, const char* newline = "\n");

		static void error(char* str, bool to_screen = true, bool to_log = true, const char* newline = "\n");
		static void warning(char* str, bool to_screen = true, bool to_log = true, const char* newline = "\n");
		static void log(char* str, bool to_screen = true, bool to_log = true, const char* newline = "\n");

		static void error(char c, bool to_screen = true, bool to_log = true, const char* newline = "\n");
		static void warning(char c, bool to_screen = true, bool to_log = true, const char* newline = "\n");
		static void log(char c, bool to_screen = true, bool to_log = true, const char* newline = "\n");
	};
}

#endif