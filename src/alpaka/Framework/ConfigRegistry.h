#ifndef Framework_ConfigRegistry_h
#define Framework_ConfigRegistry_h

#include <nlohmann/json.hpp>
#include <fstream>
#include <string>
#include <unordered_map>
#include <stdexcept>

#define INPUT_DEBUG

#ifdef INPUT_DEBUG
#include <iostream>  // for debug prints
#endif

namespace edm {

  using Config = nlohmann::json;

  class ConfigRegistry {
  public:

    ConfigRegistry() = default;

    // Just load configuration from a JSON file
    static ConfigRegistry loadFromFile(const std::string& filename) {
  #ifdef INPUT_DEBUG
      std::cout << "[ConfigRegistry] Loading configuration from file: " << filename << std::endl;
  #endif

      std::ifstream in(filename);
      if (!in.is_open()) {
        throw std::runtime_error("Cannot open configuration file: " + filename);
      }

      Config root;
      try {
        in >> root;
      } catch (const std::exception& e) {
        throw std::runtime_error(std::string("[ConfigRegistry] JSON parse failed: ") + e.what());
      }

  #ifdef INPUT_DEBUG
      std::cout << "[ConfigRegistry] Successfully parsed JSON with "
                << root.size() << " top-level entries." << std::endl;
  #endif

      return ConfigRegistry(std::move(root));
    }

    // Retrieve the config object for a given producer
const Config& getProducerConfig(const std::string& fullName) const {
#ifdef INPUT_DEBUG
  std::cout << "[ConfigRegistry] Fetching config for producer: " << fullName << std::endl;
#endif

  // Derive a normalized name: strip "alpaka_*::" if present
  std::string name = fullName;
  auto pos = name.find("alpaka_");
  if (pos != std::string::npos) {
    auto nsEnd = name.find("::", pos);
    if (nsEnd != std::string::npos) {
      name = name.substr(nsEnd + 2);
#ifdef INPUT_DEBUG
      std::cout << "[ConfigRegistry] Normalized module name: " << name << std::endl;
#endif
    }
  }

  auto it = configs_.find(name);
  if (it == configs_.end()) {
#ifdef INPUT_DEBUG
    std::cout << "[ConfigRegistry] No configuration found for '" << name
              << "', using empty config (defaults will apply)." << std::endl;
#endif
    static const Config emptyConfig = Config::object();
    return emptyConfig;
  }

#ifdef INPUT_DEBUG
  std::cout << "[ConfigRegistry] Found configuration for '" << name
            << "' with " << it->second.size() << " parameters." << std::endl;
#endif

  return it->second;
}


  private:
    explicit ConfigRegistry(Config root) {
  #ifdef INPUT_DEBUG
      std::cout << "[ConfigRegistry] Building configuration registry..." << std::endl;
  #endif

      for (auto& [key, val] : root.items()) {
  #ifdef INPUT_DEBUG
        std::cout << "  - Registering producer: " << key
                  << " (" << val.size() << " params)" << std::endl;
  #endif
        configs_.emplace(key, val);
      }

  #ifdef INPUT_DEBUG
      std::cout << "[ConfigRegistry] Registry built with " << configs_.size()
                << " producers." << std::endl;
  #endif
    }

    /// TODO: maybe it's better to have a custom class for the params
    std::unordered_map<std::string, Config> configs_;
  };

}

#endif  // Framework_ConfigRegistry_h
