import re
import argparse
import os

def eventsFromFile(file, eventMarker):
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

            if line.startswith(eventMarker):
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
    # Updates the number of particles in the event header by changing any number after ':' with numParticles
    newHeader = re.sub(r":\d+", f":{numParticles}", header)
    return newHeader


def filterParticlesPerLayer(particlesFile, hitsFile, outputFile, minLayers=4):
    genPart = eventsFromFile(particlesFile, "particles:")
    genHits = eventsFromFile(hitsFile, "hits:")

    with open(outputFile, "w", encoding="utf-8") as out:
        for event_number, (header, partEvent) in enumerate(genPart, start=1):
            try:
                _, hitEvent = next(genHits)
            except StopIteration:
                raise ValueError(
                    f"No matching hit event for particle event {event_number} "
                    f"({header}) in '{hitsFile}'."
                ) from None

            layers_by_pid = {
                line.rsplit(",", 1)[-1].strip(): set()
                for line in partEvent
            }

            for line in hitEvent:
                fields = line.split(",")
                if len(fields) != 15:
                    continue

                try:
                    layer_id = int(fields[-3])
                except ValueError:
                    continue

                particle_id = fields[-1].strip()
                if particle_id in layers_by_pid:
                    layers_by_pid[particle_id].add(layer_id)

            filtered = [
                line for line in partEvent
                if len(layers_by_pid[line.rsplit(",", 1)[-1].strip()]) >= minLayers
            ]

            out.write(updateHeader(header, len(filtered)) + "\n")
            out.writelines(line + "\n" for line in filtered)

        try:
            extra_hit_header, _ = next(genHits)
        except StopIteration:
            pass
        else:
            raise ValueError(
                f"'{hitsFile}' has extra hit events, starting with "
                f"{extra_hit_header}; the input files do not match."
            )

    print(f"File '{outputFile}' created successfully!")


# Using example
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog='extractHitsFromParquet',
                    description='Extract hits information from parquet files')

    parser.add_argument('-o', '--outFileName', default='')
    parser.add_argument('--outDir', default='.', help='Directory containing inputs and for the output file (default: current directory).')
    parser.add_argument('-p', '--isPhase2', action='store_true')
    parser.add_argument('-m', '--isMuonsOnly', action='store_true')

    args = parser.parse_args()
    outDir = args.outDir or '.'
    os.makedirs(outDir, exist_ok=True)

    phaseModifier = 'Phase1'
    if args.isPhase2: phaseModifier = 'Phase2'

    muonModifier = ''
    if args.isMuonsOnly: muonModifier = 'MuonsOnly'

    inputPartFileName = os.path.join(outDir, 'particles' + args.outFileName + muonModifier + '.txt')
    inputHitsName = os.path.join(outDir, 'hitsWithParticleID' + args.outFileName + phaseModifier + '.txt')
    outFileName = os.path.join(outDir, 'particlesFilter' + args.outFileName + muonModifier + phaseModifier + '.txt')

    filterParticlesPerLayer(inputPartFileName, inputHitsName, outFileName)

