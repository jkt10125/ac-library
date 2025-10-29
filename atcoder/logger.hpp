#ifndef ATCODER_LOGGER_HPP
#define ATCODER_LOGGER_HPP

#include "internal_type_traits.hpp"

#include <mutex>
#include <memory>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <functional>
#include <iomanip>

namespace atcoder {

struct Logger {
  public:
    static Logger &getInstance() {
        static Logger instance;
        return instance;
    }

    void setOutputFile(const std::string &logFilename, bool append = true) {
        std::lock_guard<std::mutex> globallock(myGlobalMutex);

        defaultLogFilename = logFilename;

        if (myLoggerManager.count(logFilename)) {
            return;
        }
        
        std::ios_base::openmode mode = std::ios::out;
        mode |= (append ? std::ios::app : std::ios::trunc);

        auto fileStreamPtr = std::make_unique<std::ofstream>(logFilename, mode);
        if (!fileStreamPtr->is_open()) {
            throw std::runtime_error("Logger: Unable to open log file: " + logFilename);
        }

        myLoggerManager[logFilename] = std::make_shared<LoggerData>(
            LoggerData{std::move(fileStreamPtr), std::make_shared<std::mutex>()}
        );
    }

    void log(const std::string &filename,
             const std::string &msg,
             const char *file,
             const char *func,
             const int line)
    {
        std::shared_ptr<LoggerData> loggerData;

        {   // this scoping is important for globallock to be released early
            std::lock_guard<std::mutex> globallock(myGlobalMutex);

            if (myLoggerManager.count(filename) == 0) {
                std::cerr << "Logger: Log file not set: " << filename << std::endl;
                return;
            }
            loggerData = myLoggerManager[filename];
        }

        writeToStream(loggerData, msg, file, func, line);
    }

    void log(const std::string &msg,
             const char *file,
             const char *func,
             const int line)
    {
        log(defaultLogFilename, msg, file, func, line);
    }


  private:
    struct LoggerData {
        std::unique_ptr<std::ofstream> outStreamPtr;
        std::shared_ptr<std::mutex> mutexPtr;
    };

    std::unordered_map<std::string, std::shared_ptr<LoggerData>> myLoggerManager;
    std::mutex myGlobalMutex;
    std::string defaultLogFilename;

    Logger() = default;
    ~Logger() = default;
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    void writeToStream(std::shared_ptr<LoggerData> loggerData,
                       const std::string &msg,
                       const char *file,
                       const char *func,
                       const int line)
    {
        auto now = std::chrono::system_clock::now();
        auto time = std::chrono::system_clock::to_time_t(now);
        auto micros = std::chrono::duration_cast<std::chrono::microseconds>(
                          now.time_since_epoch()) % 1'000'000;

        std::ostringstream oss;
        oss << std::put_time(std::localtime(&time), "%Y-%m-%d %H:%M:%S")
            << "." << std::setfill('0') << std::setw(6) << micros.count() << " ";

        // --- File + Line + Function ---
        oss << file << ":" << line << " (" << func << ") [ ";

        // --- Message ---
        oss << msg << " ]\n";

        std::string out = oss.str();
        std::lock_guard<std::mutex> filelock(*(loggerData->mutexPtr));
        (*(loggerData->outStreamPtr)) << out << std::flush;
    }
};

}  // namespace atcoder

// --- Macros ---
#define LOG_TO_FILE(filename, message) \
    atcoder::Logger::getInstance().log(filename, message, __FILE__, __func__, __LINE__)

#define LOG(message) \
    atcoder::Logger::getInstance().log(message, __FILE__, __func__,  __LINE__)

#endif      // ATCODER_LOGGER_HPP