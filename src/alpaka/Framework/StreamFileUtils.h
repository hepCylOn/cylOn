#ifndef Framework_StreamFileUtils_h
#define Framework_StreamFileUtils_h

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace edm::utils {

  inline std::ifstream openInputFile(const std::filesystem::path& path,
                                     std::ios::openmode mode = std::ios::binary) {

    // Expand to absolute path for clarity
    std::filesystem::path absolutePath = absolute(path);

    if (!exists(absolutePath)) {
      std::ostringstream msg;
      msg << "[FileUtils] File not found: " << absolutePath;
      throw std::runtime_error(msg.str());
    }

    std::ifstream file(absolutePath, mode);
    if (!file.is_open()) {
      std::ostringstream msg;
      msg << "[FileUtils] Failed to open file: " << absolutePath;
      throw std::runtime_error(msg.str());
    }

    file.exceptions(std::ios::failbit | std::ios::badbit);  // Enable I/O exceptions
    return file;  // RVO — no copy thanks to move semantics
  }

}  // namespace edm::utils

#endif  // Framework_StreamFileUtils_h
