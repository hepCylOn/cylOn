![alt text](assets/cylon_logo.png)

# CylOn: An(other) Geometry Agnostic Framework for Heterogenous Track Reconstruction

CylOn is a framework allowing to reconstruct HEP experiments tracks on multiple backends. It is based on the CMS Experiment: 
- [Patatrack](https://patatrack.web.cern.ch/patatrack/wiki/) pixel tracking CA-based chain;
- Alpaka-based heterogenous framework, with a simplified version. 

> [!IMPORTANT] 
> The vast majority of the algorithmic backbone comes from the work done in the original [Patatrack standalone pixel tracking repository](https://github.com/hpc4lhc/pixeltrack-standalone), from which this repository is forked, and the furter improvements brought to the framework in [cms-sw](https://github.com/cms-sw/cmssw). The full merit and credit goes to the developers there (some of which are also working on this project). Find here a list of references
> 
> 

## Setup

From the framework's root, use the command `make alpaka` to build the framework. Option `-j` can be used to parallelize the build. Check https://github.com/hepCylOn/cylOn/tree/checkpointFromHits#troubleshooting for more information on the needed definitions and possible issues that can happen at this step.

After setting up, set the environment with `. env.sh`. Then, create a virtual python environment to download the needed python libraries:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

After the intial setup, to use the framework just use the following commands from the `cylOn` directory:

```sh
. env.sh
source .venv/bin/activate
```

To deactivate the python environment use the command `deactivate`.

## Running with `cylon.sh`

Run the script from the `cylOn` project directory. It builds the `alpaka` target with `make` and then runs `./alpaka` with the selected options. By default it uses the CUDA backend, processes up to 100 events, and builds with 8 jobs.

```bash
./cylon.sh \
	--backend cuda \
	--max-events 10 \
	--from-hits \
	--phase2 \
	--collider-ml \
	--tracker-type PixelOnly \
	--validation
```

`--backend` accepts `cuda`, `serial`, or `rocm`. The reconstruction options map to the corresponding `alpaka` flags: `--from-hits` to `--fromHits`, `--phase2` to `--isPhase2`, `--collider-ml` to `--isColliderML`, and `--validation` to `--validation`. `--tracker-type` has to be used with `--collider-ml` and accepts `PixelOnly`, `PixelPlusShortStrips`, or `AllTracker`.

Other useful options include:

- `--max-events N` sets the event limit.
- `--make-jobs N` sets the parallelism used by `make`.
- `--run-twice` runs the reconstruction twice; `--first-extra ARG` and `--second-extra ARG` add arguments to the respective runs, and can be included multiple times for distinct arguments.
- `--autograph --input FILE` runs `scripts/autograph.py` before the build. Use `--python-script SCRIPT` to select another script and `--python-arg ARG` to pass it extra arguments. Can be included multiple times for multiple arguments. `--input` is for the Python script and is not passed to `alpaka`.

For example, this wrapper invocation builds and runs the equivalent of `./alpaka --cuda --maxEvents 10 --fromHits --isPhase2 --isColliderML --trackerType PixelOnly --validation`:

```bash
./cylon.sh --backend cuda --max-events 10 --from-hits --phase2 --collider-ml --tracker-type PixelOnly --validation
```

### Running twice

Add `--run-twice` to run the same reconstruction twice with the common options. The script compiles `alpaka` once before both runs:

```bash
./cylon.sh --backend cuda --max-events 10 --from-hits --phase2 --collider-ml --tracker-type PixelOnly --validation --run-twice
```

To vary an option between runs, pass each token of that option through the corresponding extra flag. For example, to run once with one thread and once with `--runSimTracks` and once without:

```bash
./cylon.sh \
	--backend cuda \
	--max-events 10 \
	--from-hits \
	--phase2 \
	--collider-ml \
	--tracker-type PixelOnly \
	--validation \
	--run-twice \
	--first-extra --runSimTracks
```

This runs the equivalent of:

```bash
./alpaka --cuda --maxEvents 10 --fromHits --isPhase2 --isColliderML --validation --trackerType PixelOnly --runSimTracks
./alpaka --cuda --maxEvents 10 --fromHits --isPhase2 --isColliderML --validation --trackerType PixelOnly
```

For completeness, an example with parsing an argument that needs a value would be:

```bash
./cylon.sh \
	--backend cuda \
	--from-hits \
	--phase2 \
	--collider-ml \
	--tracker-type PixelOnly \
	--validation \
	--run-twice \
	--first-extra --numberOfThreads \
	--first-extra 1 \
	--second-extra --numberOfThreads \
	--second-extra 4
```

This runs the equivalent of:

```bash
./alpaka --cuda --maxEvents 10 --fromHits --isPhase2 --isColliderML --validation --trackerType PixelOnly --numberOfThreads 1
./alpaka --cuda --maxEvents 10 --fromHits --isPhase2 --isColliderML --validation --trackerType PixelOnly --numberOfThreads 4
```

## Producing input data from ColliderML

-----
<details>

<summary>Citations for the ColliderML dataset used</summary>


https://arxiv.org/abs/2512.15230

```
@misc{elitez2025collidermlreleaseopendatadetectorhighluminosity,
      title={ColliderML: The First Release of an OpenDataDetector High-Luminosity Physics Benchmark Dataset}, 
      author={Doğa Elitez and Paul Gessinger and Daniel Murnane and Marcus Selchou Raaholt and Andreas Salzburger and Stine Kofoed Skov and Andreas Stefl and Anna Zaborowska},
      year={2025},
      eprint={2512.15230},
      archivePrefix={arXiv},
      primaryClass={hep-ex},
      url={https://arxiv.org/abs/2512.15230}, 
}
```

```
@dataset{colliderml_release1_2025,
  title={{ColliderML Dataset Release 1}},
  author={{ColliderML Collaboration}},
  year={2025},
  publisher={Hugging Face},
  howpublished={\url{https://huggingface.co/datasets/OpenDataDetector/ColliderML-Release-1}},
  note={Simulation performed using ACTS and the Open Data Detector}
}
```

</details>

-----

All the scripts used to generate ColliderML data can be found in the `scripts` folder.

### Producing hits information and map from hits to particles association

Script [extractHitsFromParquet.py](./scripts/extractHitsFromParquet.py) can be used to extract the hits information from the ColliderML dataset. It will save the hits in a .txt file that will contain the following info:

```
events:NUMBER_OF_EVENTS
hits:NUMBER_OF_HITS_EV1
modules:NUMBER_OF_LAYERS
hit1_ph,hit1_ph,hit1_dummy_error,hit1_dummy_error,hit1_xg,hit1_yg,hit1_zg,hit1_rg,hit1_phiShort,hit1_ph,hit1_ph,hit1_ph,hit1_layerID,hit1_partID
...
hitN_ph,hitN_ph,hitN_dummy_error,hitN_dummy_error,hitN_xg,hitN_yg,hitN_zg,hitN_rg,hitN_phiShort,hitN_ph,hitN_ph,hitN_ph,hitN_layerID,hitN_partID
0,nHits_layer(0),nHits_layer(1)+prev,...,nHitsLayer(N-1)+prev,nHitsLayer(N)+prev
hits:NUMBER_OF_HITS_EV2
...
```

The output .txt file can have a large size, depending on the amount of events size. It can be compressed with command: `gzip -c hits.txt > hits.txt.gz`, and `hits.txt.gz` can be directly used as an input to CylOn. The `-c` option keeps the input file.

The map from hits to particles is also saved as a .txt file with the number of events as `events:NUMBER_OF_EVENTS`, the number of hits per event as `hits:NUMBER_OF_HITS`, and the partID associated to each hit, in the same order that they appear in the hits file.

### Producing particles information

Script [extractParticlesFromParquet.py](./scripts/extractParticlesFromParquet.py) can be used to extract the particles information from the ColliderML dataset. It will save the particles in a .txt file that will contain the following info:

```
particles:NUMBER_OF_PARTICLES_EV1
part1_vx,part1_vy,part1_vz,part1_px,part1_py,part1_pz,part1_e,part1_pt,part1_eta,part1_phi,part1_m,part1_c,part1_pid,part1_id
...
partN_vx,partN_vy,partN_vz,partN_px,partN_py,partN_pz,partN_e,partN_pt,partN_eta,partN_phi,partN_m,partN_c,partN_pid,partN_id
particles:NUMBER_OF_PARTICLES_EV2
...
```
There is a pre-selection applied to the particles at this stage in terms of the particles pT (> 0.9 GeV), |eta| (< 3.2), |dz| (< 30 cm), dxy (< 2.5 cm) and |pid| (optional, set to the muon PDGID = 13).

Afterwards, script [extractParticlesWithMoreHits.py](./scripts/extractParticlesWithMoreHits.py) can be used to only keep the particles associated to at least 4 hits, giving a better chance to the reconstruction algorithm to reconstruct these particles.

-----

<details>

<summary>Example of how to run the scripts</summary>

### Simple test workflow

Run the commands below to produce some testing files:

```sh
# From the cylon directory do:

python3 scripts/extractHitsFromParquet.py -o AllTracker --outDir data -p -s -l -n 10
python3 scripts/extractParticlesFromParquet.py -o AllTracker --outDir data -n 10
python3 scripts/extractParticlesWithMoreHits.py -o AllTracker --outDir data -p

# Don't use the -c option below if you want to delete the input files

gzip -c data/hitsWithParticleIDAllTrackerPhase2.txt > data/hitsWithParticleIDAllTrackerPhase2.txt.gz
gzip -c data/mapHitsToParticlesAllTrackerPhase2.txt > data/mapHitsToParticlesAllTrackerPhase2.txt.gz
gzip -c data/particlesFilterAllTrackerPhase2.txt > data/particlesFilterAllTrackerPhase2.txt.gz
```

The configs should be already set so that `./cylon.sh --backend cuda --max-events 10 --from-hits --phase2 --collider-ml --tracker-type AllTracker --validation` can be executed. If not, just change [testColliderMLAllTrackerPhase2FromHits.json](./configs/testColliderMLAllTrackerPhase2FromHits.json) to the correct file names.

</details>

-----

## Updating a config's CA parameters

The [autograph.py](./scripts/autograph.py) script extracts some of the CA parameters from the geometry defined by the input hits true position in the detector. To apply the newly extracted parameters, that are saved in a file named `line_sweep_edges.txt`, use the script [modifyConfigsCAParams.py](./scripts/modifyConfigsCAParams.py) as shown below to create a new config file with these parameters:

```bash
python3 scripts/modifyConfigsCAParams.py -c configs/templateConfig.json -i line_sweep_edges.txt -o configs/OUTPUT_FILE_NAME.json
```

Running with `SimPixelTracks` will produce plots that are used to extract updated CA parameters by the [plot-simPixelTracks.py](./scripts/plot-simPixelTracks.py) script in a .txt file. Apply them to a new config based on the output config of the step above:

```bash
python3 scripts/modifyConfigsCAParams.py -c configs/OUTPUT_FILE_NAME.json -i SIMPIXELTRACKS_OUTPUT.txt -o configs/OUTPUT_FILE_NAME.json -f -s
```

The `-s` option tells the script that it is using the output of the SimPixelTracks, which will modify distinct modules and entries of the config json.

## Troubleshooting

Executing `./cylon.sh` also executes the `Makefile` which can show some issues. Here are some of the problems that can appear (non-exhaustive list)

### Proxy for curl

The framework uses curl to download some files that are needed. Sometimes the proxy `-x` option is needed with the correct proxy identification. The proxy can be defined either as:

```sh
export MYPROXY=http://proxy.example:8080
```

Or parsing it when running the make commmand:

```sh
make -j 8 MYPROXY=http://proxy.example:8080 alpaka
```

If no proxy is needed, leave `MYPROXY` unchanged.

### Checking the GCC toolchain

If the compiler cannot find `cmath` it will be probably related to the `GCC_TOOLCHAIN` variable defined in the `Makefile`.

`GCC_TOOLCHAIN` is a Make variable, so checking it with `echo "$GCC_TOOLCHAIN"` in your shell may show an empty value. To see the value Make will use, run this from the project directory:

```sh
make --eval='print-toolchain:; @echo "GCC_TOOLCHAIN=$(GCC_TOOLCHAIN)"' print-toolchain
```

This only parses the Makefile and runs the temporary `print-toolchain` target; it does not build `alpaka`.

Check that the printed path is the GCC installation root expected by your compiler. If automatic detection picks the wrong path—for example, because CXX points to a ccache wrapper—you can override it for a build:

```sh
make -j 8 GCC_TOOLCHAIN=/path/to/gcc/toolchain alpaka
```

Or define it before building as:

```sh
export GCC_TOOLCHAIN=/path/to/gcc/toolchain
```
