# -*- coding: utf-8 -*-
"""
Filters particles from the particles.txt file, keeping only
those that have hits in at least 4 different layers,
processing event by event to save memory.
The output file maintains the same structure as the input file,
and updates the number of particles in the header of each event.
"""

import re
import argparse

def eventsFromFile(file, eventMarker):
    """
    Reads a text file and yields events one by one.

    Returns:
        (header, eventLines)
    """

    currentEvent = []
    header = None

    with open(file, "r", encoding="utf-8") as f:

        for rawLine in f:

            line = rawLine.strip()

            if not line:
                continue

            # Ignore global metadata
            if line.startswith("events:"):
                continue

            # New event found
            if line.startswith(eventMarker):

                # Yield previous event BEFORE starting new one
                if header is not None:
                    yield header, currentEvent

                header = line
                currentEvent = []

            else:
                currentEvent.append(line)

        # Yield final event
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


def filterParticlesPerLayer(particlesFile, hitsFile, outputFile, minLayers=4):
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

            EXPECTED_HIT_COLUMNS = 15  # adjust for the hits format

            for line in hitEvent:

                pieces = line.strip().split(',')

                # Ignore summary lines or malformed lines
                if len(pieces) != EXPECTED_HIT_COLUMNS:
                    continue
                
                try:
                    layerId = int(pieces[-3])
                    particleID = pieces[-1].strip()

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
    parser = argparse.ArgumentParser(
                    prog='extractHitsFromParquet',
                    description='Extract hits information from parquet files')

    parser.add_argument('-o', '--outFileName', default='')
    parser.add_argument('-p', '--isPhase2', action='store_true')
    parser.add_argument('-m', '--isMuonsOnly', action='store_true')

    args = parser.parse_args()

    phaseModifier = 'Phase1'
    if args.isPhase2: phaseModifier = 'Phase2'

    muonModifier = ''
    if args.isMuonsOnly: muonModifier = 'MuonsOnly'

    inputPartFileName = 'particles' + args.outFileName + muonModifier + '.txt'
    inputHitsName = 'hitsWithParticleID' + args.outFileName + phaseModifier + '.txt'
    outFileName = 'particlesFilter' + args.outFileName + muonModifier + phaseModifier + '.txt'

    filterParticlesPerLayer(inputPartFileName, inputHitsName, outFileName)


