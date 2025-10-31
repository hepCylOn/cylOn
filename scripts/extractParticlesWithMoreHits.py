# -*- coding: utf-8 -*-
"""
Filters particles from the particles.txt file, keeping only
those that have hits in at least 4 different layers,
processing event by event to save memory.
The output file maintains the same structure as the input file,
and updates the number of particles in the header of each event.
"""

import re

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
                if currentEvent:
                    yield header, currentEvent
                    currentEvent = []
                header = line
            else:
                currentEvent.append(line)
        if currentEvent:
            yield header, currentEvent


def updateHeader(header, numParticles):
    """
    Updates the number of particles in the event header. 
    Example: 'particles:10' -> 'particles:5'
    """
    # Substitutes any number after ':' by numParticles
    newHeader = re.sub(r":\d+", f":{numParticles}", header)
    return newHeader


def filterParticlesPerLayer(particlesFile, hitsFile, outputFile="particlesFilter.txt", minLayers=4):
    """
    Filters particles, keeping only those that have hits in at least 'minLayers' layers. 
    It maintains the same structure as the input file, including empty events (header only),
    and updates the header of each event with the number of filtered particles.
    """
    genPart = eventsFromFile(particlesFile, "particles:")
    genHits = eventsFromFile(hitsFile, "hits:")

    with open(outputFile, "w", encoding="utf-8") as out:
        for header, partEvent in genPart:
            _, hitEvent = next(genHits)

            # Maps particleID -> set of layerIDs
            particleIDs = [line.split(",")[-1].strip() for line in partEvent]
            layersPerParticle = {pid: set() for pid in particleIDs}

            # Variable to skip the line after 'module:'
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
                    particleID = pieces[-1]
                except ValueError:
                    continue
                if particleID in layersPerParticle:
                    layersPerParticle[particleID].add(layerId)

            # Selects particles that pass the criteria (at least 4 hits in distinct layers)
            filteredParticles = []
            for line in partEvent:
                pid = line.split(",")[-1].strip()
                layers = sorted(layersPerParticle.get(pid, set()))
                if len(layers) >= minLayers:
                    modifiedLine = line
                    filteredParticles.append(modifiedLine)

            # Updates the header with the new number of particles
            newHeader = updateHeader(header, len(filteredParticles))
            out.write(newHeader + "\n")

            for line in filteredParticles:
                out.write(line + "\n")

    print(f"File '{outputFile}' created successfully!")


# Using example
if __name__ == "__main__":
    filterParticlesPerLayer("/data/user/borzari/cmssw/pixeltrack-standalone/data/particles.txt", "/data/user/borzari/cmssw/pixeltrack-standalone/data/hits.txt")


