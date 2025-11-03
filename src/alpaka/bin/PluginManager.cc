#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>

#include "PluginManager.h"

#ifndef LIB_DIR
#error "LIB_DIR undefined"
#endif

#define STR_EXPAND(x) #x
#define STR(x) STR_EXPAND(x)

#define FW_DEBUG

namespace edmplugin {
  PluginManager::PluginManager() {
    std::ifstream pluginMap(STR(LIB_DIR) "/plugins.txt");
    std::string plugin, library;
    while (pluginMap >> plugin >> library) {
#ifdef FW_DEBUG
      std::cout << "plugin " << plugin << " in " << library << std::endl;
      /// TODO: add error if plugin has no library
#endif
      pluginToLibrary_[plugin] = library;
    }
  }

  SharedLibrary const& PluginManager::load(std::string const& pluginName) {
    std::lock_guard<std::recursive_mutex> guard(mutex_);

    auto libName = pluginToLibrary_.at(pluginName);

    auto found = loadedPlugins_.find(libName);
    if (found == loadedPlugins_.end()) {
      auto ptr = std::make_shared<SharedLibrary>(STR(LIB_DIR) "/" + libName);
      loadedPlugins_[libName] = ptr;
      return *ptr;
    }
    return *(found->second);
  }
}  // namespace edmplugin
