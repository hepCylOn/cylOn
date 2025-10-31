# -*- coding: utf-8 -*-
"""
Filters particles from the particles.txt file, keeping only
those that have hits in at least 4 different layers,
processing event by event to save memory.
The output file maintains the same structure as the input file,
and updates the number of particles in the header of each event.
"""

import re
import itertools
from itertools import zip_longest


def eventsFromFile(file, eventMarker):
    """
    A generator that reads a text file and returns events one by one.
    Each event is a tuple (header, list of particle/hit lines)
    """
    currentEvent = []
    header = None

    with open(file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(eventMarker):
                # Return the event even if empty
                if header is not None:
                    yield header, currentEvent
                    currentEvent = []
                header = line
            else:
                currentEvent.append(line)

        # At the end of the file returns the event, even if empty
        if header is not None:
            yield header, currentEvent


def updateHeader(header, numParticles):
    """
    Updates the number of particles in the event header. 
    Example: 'particles:10' -> 'particles:5'
    """
    # Substitutes any number after ':' by numParticles
    newHeader = re.sub(r":\d+", f":{numParticles}", header)
    return newHeader


def filterParticlesPerLayer(
    particlesFile,
    hitsFile,
    outputParticlesFile="particlesLastWithMinLayers.txt",
    outputHitsFile="hitsLastWithMinLayers.txt",
    minLayers=4,
    numLayers=18
):
    """
    Keeps only the last particle of each event, provided it has hits
    in at least `minLayers` distinct layers (0-17).
    Generates synchronized particles and hits files with updated headers.
    """
    genPart = eventsFromFile(particlesFile, "particles:")
    genHits = eventsFromFile(hitsFile, "hits:")

    with open(outputParticlesFile, "w", encoding="utf-8") as outPart, \
         open(outputHitsFile, "w", encoding="utf-8") as outHits:

        for (header, partEvent), (hitHeader, hitEvent) in zip_longest(genPart, genHits, fillvalue=(None, [])):
            if header is None:
                break

            # Event without particles
            if len(partEvent) == 0:
                outPart.write(updateHeader(header, 0) + "\n")
                outHits.write(f"hits:0\nmodule:{numLayers}\n" + ",".join(["0"] * (numLayers + 1)) + "\n")
                continue

            # Gets only latest muon in event (highest pT)
            lastParticleLine = partEvent[-1]
            lastParticleID = lastParticleLine.split(",")[-1].strip()

            # Maps layers number with hits
            layersForLastParticle = set()

            skipNextLine = False
            for line in hitEvent:
                if skipNextLine:
                    skipNextLine = False
                    continue
                if line.startswith("module:"):
                    skipNextLine = True
                    continue

                pieces = line.strip().split(',')
                if len(pieces) < 2:
                    continue
                try:
                    layerId = int(pieces[-2])
                    particleID = pieces[-1].strip()
                except ValueError:
                    continue

                if particleID == lastParticleID:
                    layersForLastParticle.add(layerId)

            # Check if at least 4 layers
            if len(layersForLastParticle) < minLayers:
                # If particle doesn't match the requisite, still writes an empty event
                outPart.write(updateHeader(header, 0) + "\n")
                outHits.write(f"hits:0\nmodule:{numLayers}\n" + ",".join(["0"] * (numLayers + 1)) + "\n")
                continue

            # If it pass, write the particle
            outPart.write(updateHeader(header, 1) + "\n")
            outPart.write(lastParticleLine + "\n")

            # Filter matching hits
            filteredHits = []
            skipNextLine = False
            for line in hitEvent:
                if skipNextLine:
                    skipNextLine = False
                    continue
                if line.startswith("module:"):
                    skipNextLine = True
                    continue

                pieces = line.strip().split(',')
                if len(pieces) < 2:
                    continue
                particleID = pieces[-1].strip()
                if particleID == lastParticleID:
                    filteredHits.append(line)

            # Calculates number of hits per layer
            layerCounts = [0] * numLayers
            for line in filteredHits:
                pieces = line.strip().split(',')
                try:
                    layerId = int(pieces[-2])
                except ValueError:
                    continue
                if 0 <= layerId < numLayers:
                    layerCounts[layerId] += 1

            # Cumulative calculation of module hits
            cumulative = [0] + list(itertools.accumulate(layerCounts))

            # Writes filtered hits with updated modules
            outHits.write(f"hits:{len(filteredHits)}\n")
            for line in filteredHits:
                outHits.write(line + "\n")
            outHits.write(f"module:{numLayers}\n")
            outHits.write(",".join(map(str, cumulative)) + "\n")

    print(f"Successfully generated files:")
    print(f"  - {outputParticlesFile}")
    print(f"  - {outputHitsFile}")


# Using example
if __name__ == "__main__":
    filterParticlesPerLayer("particlesMuonsOnly.txt", "hits.txt")
