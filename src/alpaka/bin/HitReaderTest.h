#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdint>

#include "AlpakaDataFormats/TrackingRecHitsHost.h"

namespace hitReaderTest {

    struct Event {
        uint32_t nHits;
        uint32_t nModules;
    };

    inline void skip_separator(const char*& p) {
        while (*p == ' ' || *p == ',' || *p == '\t')
            ++p;
    }

    inline float parse_float(const char*& p) {
        char* end;
        float value = strtof(p, &end);
        p = end;
        skip_separator(p);
        return value;
    }

    inline int parse_int(const char*& p) {
        char* end;
        int value = strtol(p, &end, 10);
        p = end;
        skip_separator(p);
        return value;
    }

    inline unsigned parse_uint(const char*& p) {
        char* end;
        unsigned value = strtoul(p, &end, 10);
        p = end;
        skip_separator(p);
        return value;
    }

    inline reco::TrackingRecHitHost readEvent(std::ifstream& file) {

        Event event;

        reco::TrackingRecHitHost auxRecHitHost(cms::alpakatools::host(), 0, 0);

        std::string line;

        // Skip empty lines
        while (std::getline(file, line)) {
            if (!line.empty())
                break;
        }

        if (file.eof())
            return auxRecHitHost;

        if (line.rfind("hits:", 0) != 0) {
            std::cerr << "Expected hits header\n";
            return auxRecHitHost;
        }

        event.nHits = std::strtoul(line.c_str() + 5, nullptr, 10);

        if (!std::getline(file, line)) {
            std::cerr << "Missing module header\n";
            return auxRecHitHost;
        }

        if (line.rfind("module:", 0) != 0) {
            std::cerr << "Expected module header\n";
            return auxRecHitHost;
        }

        event.nModules = std::strtoul(line.c_str() + 7, nullptr, 10);

        reco::TrackingRecHitHost recHitHost(cms::alpakatools::host(), event.nHits, event.nModules);

        auto hitsView = recHitHost.view<reco::TrackingRecHitSoA>();
        auto modsView = recHitHost.view<reco::HitModuleSoA>();

        for (uint32_t i = 0; i < event.nHits; ++i) {

            if (!std::getline(file, line)) {
                std::cerr << "Unexpected EOF while reading hits\n";
                return auxRecHitHost;
            }

            const char* p = line.c_str();

            hitsView[i].xLocal() = parse_float(p);
            hitsView[i].yLocal() = parse_float(p);
            hitsView[i].xerrLocal() = parse_float(p);
            hitsView[i].yerrLocal() = parse_float(p);

            hitsView[i].xGlobal() = parse_float(p);
            hitsView[i].yGlobal() = parse_float(p);
            hitsView[i].zGlobal() = parse_float(p);
            hitsView[i].rGlobal() = parse_float(p);

            hitsView[i].iphi() = static_cast<int16_t>(parse_int(p));

            hitsView[i].chargeAndStatus().charge = parse_float(p);

            hitsView[i].clusterSizeX() = static_cast<int16_t>(parse_int(p));
            hitsView[i].clusterSizeY() = static_cast<int16_t>(parse_int(p));

            hitsView[i].detectorIndex() = static_cast<uint16_t>(parse_uint(p));

            hitsView.offsetBPIX2() = static_cast<int32_t>(parse_int(p));
        }

        if (!std::getline(file, line)) {
            std::cerr << "Missing moduleStart line\n";
            return auxRecHitHost;
        }

        {
            const char* p = line.c_str();

            for (uint32_t i = 0; i < event.nModules + 1; ++i) {
                modsView[i].moduleStart() = parse_uint(p);
            }
        }

        return recHitHost;
    }

}

namespace mapReaderTest {

    struct Event {
        uint32_t nHits;
    };

    inline unsigned parse_uint(const char*& p) {
        char* end;
        unsigned value = strtoul(p, &end, 10);
        p = end;
        return value;
    }

    inline utils::SimpleMapHost readEvent(std::ifstream& file) {

        Event event;

        utils::SimpleMapHost auxMapHost(0, cms::alpakatools::host());

        std::string line;

        // Skip empty lines
        while (std::getline(file, line)) {
            if (!line.empty())
                break;
        }

        if (file.eof())
            return auxMapHost;

        if (line.rfind("hits:", 0) != 0) {
            std::cerr << "Expected hits header\n";
            return auxMapHost;
        }

        event.nHits = std::strtoul(line.c_str() + 5, nullptr, 10);

        utils::SimpleMapHost mapHost(event.nHits, cms::alpakatools::host());

        auto mapView = mapHost.view();

        for (uint32_t i = 0; i < event.nHits; ++i) {

            if (!std::getline(file, line)) {
                std::cerr << "Unexpected EOF while reading hits\n";
                return auxMapHost;
            }

            const char* p = line.c_str();

            mapView[i].id() = parse_uint(p);

        }

        return mapHost;
    }

}