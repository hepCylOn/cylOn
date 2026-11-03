![alt text](assets/cylon_logo.png)

# CylOn: An(other) Geometry Agnostic Framework for Heterogenous Track Reconstruction

CylOn is a framework allowing to reconstruct HEP experiments tracks on multiple backends. It is based on the CMS Experiment: 
- [Patatrack](https://patatrack.web.cern.ch/patatrack/wiki/) pixel tracking CA-based chain;
- Alpaka-based heterogenous framework, with a simplified version. 

> [!IMPORTANT] 
> The vast majority of the algorithmic backbone comes from the work done in the original [Patatrack standalone pixel tracking repository](https://github.com/hpc4lhc/pixeltrack-standalone), from which this repository is forked, and the furter improvemnts brought to the framework in [cms-sw](https://github.com/cms-sw/cmssw). The full merit and credit goes to the developers there (some of which are also working on this project). Find here a list of references
> 
> 

## Setup

Use the command `make alpaka` to build the framework. Option `-j` can be used to parallelize the build. Check https://github.com/hepCylOn/cylOn/tree/checkpointFromHits#troubleshooting for more information on the needed definitions and possible issues that can happen at this step.

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

This is useful, since running with `SimPixelTracks` extracts parameters for the CA algorithm. It still needs to be parsed to the config in an automated way.

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
